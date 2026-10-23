"""
Keep the g-factor data in step between the Zenodo draft (where DataPipe saves new sessions), the Sussex server
(gfactor/tasks) and a local folder. Safe to run as often as you like, e.g. daily from launchd.

Each run, under --local (default ~/gfactor):
  1. zenodo_raw/      mirror of the Zenodo draft's files (loose .json and DataPipe's zips); only new or changed
                      files are downloaded
  2. zenodo_checked/  those files after the study's own PHP checks (import_datapipe.py), in the Sussex format
  3. sussex_tasks/    a fresh copy of gfactor/tasks from the server (one tar stream over SSH)
  4. checked files the server doesn't have yet are copied to gfactor/tasks. Nothing on the server is ever
     overwritten: a file whose name is already there with different contents goes to conflicts/ instead.
  5. gfactor_data/    everything together - the complete dataset. For credits:
                      credit_report_gfactor.py --no-download --json-dir ~/gfactor/gfactor_data
If the server can't be reached, steps 3-4 are skipped and gfactor_data/ still gets the new files (exit code 2).

Usage:
  python gfactor_sync.py              # full sync
  python gfactor_sync.py --dry-run    # everything except writing to the Sussex server

Settings come from the environment or ~/.gfactor_sync.env (KEY=value lines; keep it chmod 600):
  ZENODO_TOKEN     Zenodo personal access token with only the deposit:write scope (reads drafts; can't publish)
  ZENODO_DRAFT_ID  the draft's record ID, from its zenodo.org address
Needs `requests` (e.g. PhD-relaunch's .venv), PHP on the PATH, and SSH key access to unix.sussex.ac.uk.
Prints counts only, never participant data.
"""
import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import tarfile
from datetime import datetime
from pathlib import Path

import requests

from import_datapipe import SSH_HOST, REMOTE_DIR, check_files

ZENODO_API = os.environ.get('ZENODO_API', 'https://zenodo.org/api')  # overridable for testing
ENV_FILE = Path.home() / '.gfactor_sync.env'


def load_settings():
    """Environment variables win; otherwise KEY=value lines from ~/.gfactor_sync.env."""
    settings = {}
    if ENV_FILE.exists():
        for line in ENV_FILE.read_text().splitlines():
            if '=' in line and not line.lstrip().startswith('#'):
                key, value = line.split('=', 1)
                settings[key.strip()] = value.strip().strip('"\'')
    settings.update({k: v for k, v in os.environ.items() if k in ('ZENODO_TOKEN', 'ZENODO_DRAFT_ID')})
    missing = [k for k in ('ZENODO_TOKEN', 'ZENODO_DRAFT_ID') if not settings.get(k)]
    if missing:
        sys.exit(f'Missing {", ".join(missing)} - set them in the environment or {ENV_FILE}')
    return settings


def md5(path):
    digest = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def mirror_zenodo(token, draft_id, raw_dir):
    """Download the draft's files that are new or changed since the last run. Returns (downloaded, unchanged)."""
    session = requests.Session()
    session.headers['Authorization'] = f'Bearer {token}'
    response = session.get(f'{ZENODO_API}/deposit/depositions/{draft_id}', timeout=60)
    if response.status_code in (401, 403):
        sys.exit(f'Zenodo refused the token ({response.status_code}) - check ZENODO_TOKEN and its deposit:write scope')
    response.raise_for_status()
    deposition = response.json()
    if deposition.get('submitted') or deposition.get('state') == 'done':
        print('Warning: this Zenodo record has been published - raw participant data may be public')
    downloaded = unchanged = 0
    for entry in deposition.get('files', []):
        name = Path(entry.get('filename') or entry.get('key')).name  # Zenodo records have no folders
        expected = (entry.get('checksum') or '').replace('md5:', '')
        local = raw_dir / name
        if local.exists() and expected and md5(local) == expected:
            unchanged += 1
            continue
        url = (entry.get('links') or {}).get('download') or f'{deposition["links"]["bucket"]}/{name}'
        partial = local.with_name(local.name + '.part')
        with session.get(url, stream=True, timeout=300) as download:
            download.raise_for_status()
            with open(partial, 'wb') as f:
                for chunk in download.iter_content(1 << 20):
                    f.write(chunk)
        if expected and md5(partial) != expected:
            partial.unlink()
            sys.exit(f'Checksum mismatch downloading {name} from Zenodo - try again later')
        partial.replace(local)
        downloaded += 1
    return downloaded, unchanged


def pull_sussex(tasks_dir):
    """Replace tasks_dir with a fresh copy of gfactor/tasks, streamed as one tar.gz over SSH. Returns file count."""
    fresh = tasks_dir.with_name(tasks_dir.name + '.new')
    shutil.rmtree(fresh, ignore_errors=True)
    fresh.mkdir(parents=True)
    parent, folder = REMOTE_DIR.rstrip('/').rsplit('/', 1)
    stream = subprocess.Popen(['ssh', '-o', 'BatchMode=yes', SSH_HOST, f"tar czf - -C '{parent}' '{folder}'"],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    count = 0
    try:
        with tarfile.open(fileobj=stream.stdout, mode='r|gz') as tar:
            for member in tar:
                name = Path(member.name).name  # only the file's own name, so nothing lands outside the folder
                if member.isfile() and name.endswith('.json'):
                    (fresh / name).write_bytes(tar.extractfile(member).read())
                    count += 1
    except tarfile.ReadError:
        pass  # empty or broken stream - reported below with ssh's own error
    status = stream.wait()
    if status not in (0, 1) or not count:  # GNU tar exits 1 when a file changed while it was being read
        raise RuntimeError(f'could not read {REMOTE_DIR} (ssh/tar exit {status}): {stream.stderr.read().decode().strip()[:200]}')
    shutil.rmtree(tasks_dir, ignore_errors=True)
    fresh.rename(tasks_dir)
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--local', type=Path, default=Path.home() / 'gfactor', help='local folder (default ~/gfactor)')
    parser.add_argument('--dry-run', action='store_true', help="don't write anything to the Sussex server")
    args = parser.parse_args()
    print(f'--- g-factor sync {datetime.now():%Y-%m-%d %H:%M}' + (' (dry run)' if args.dry_run else ''))
    settings = load_settings()
    if not shutil.which('php'):
        sys.exit("PHP isn't on the PATH - it's needed to run the study's own checks")

    root = args.local.expanduser().resolve()
    raw, checked_dir, tasks, data, conflicts = (root / d for d in
                                                ('zenodo_raw', 'zenodo_checked', 'sussex_tasks', 'gfactor_data', 'conflicts'))
    for d in (raw, checked_dir, data, conflicts):
        d.mkdir(parents=True, exist_ok=True)

    # 1-2. Zenodo -> local, then the study's own checks
    downloaded, unchanged = mirror_zenodo(settings['ZENODO_TOKEN'], settings['ZENODO_DRAFT_ID'], raw)
    print(f'Zenodo: {downloaded} file(s) downloaded, {unchanged} unchanged')
    checked, rejected, ignored = check_files(raw, checked_dir)
    print(f'Checked: {len(checked)} passed, {len(rejected)} rejected, {len(ignored)} ignored (not study files)')
    for name, reason in rejected:
        print(f'  rejected {name}: {reason}')  # study file names are participant IDs, not personal details

    # 3-4. Sussex <-> local
    server_ok = True
    try:
        on_server = pull_sussex(tasks)
        print(f'Sussex: {on_server} file(s) in {REMOTE_DIR}')
        to_upload, conflicted = [], 0
        for path in checked:
            existing = tasks / path.name
            if not existing.exists():
                to_upload.append(path)
            elif existing.read_bytes() != path.read_bytes():
                shutil.copy2(path, conflicts / path.name)
                conflicted += 1
        if conflicted:
            print(f'  {conflicted} file(s) differ from the server copy - left alone, Zenodo versions in {conflicts}')
        if to_upload and args.dry_run:
            print(f'  dry run: would copy {len(to_upload)} new file(s) to the server')
        elif to_upload:
            subprocess.run(['scp', '-q', '-o', 'BatchMode=yes', *map(str, to_upload), f'{SSH_HOST}:{REMOTE_DIR}/'],
                           check=True, capture_output=True)
            for path in to_upload:
                shutil.copy2(path, tasks / path.name)  # the local copy of the server now matches it
            print(f'  copied {len(to_upload)} new file(s) to the server')
        else:
            print('  nothing new to copy to the server')
    except (RuntimeError, subprocess.CalledProcessError, tarfile.TarError, OSError) as e:
        server_ok = False
        print(f'Warning: Sussex step skipped - {getattr(e, "stderr", None) or e}'.strip())

    # 5. Everything together: the server's files, plus checked Zenodo files the server doesn't have
    sources = sorted(tasks.glob('*.json')) if tasks.exists() else []
    sources += [path for path in checked if not (tasks / path.name).exists()]
    for path in sources:
        shutil.copy2(path, data / path.name)
    print(f'Local: {len(list(data.glob("*.json")))} file(s) in {data}')
    sys.exit(0 if server_ok else 2)


if __name__ == '__main__':
    main()
