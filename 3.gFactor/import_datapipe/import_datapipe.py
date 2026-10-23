"""
Check g-factor study files collected through DataPipe, using the study's own PHP checks.

While the study runs from GitHub Pages, each task saves to DataPipe (see methods/save.js), which stores the files in
Zenodo or Google Drive. This script runs every file through that task's original PHP script via run_php.php, so the
checked files are exactly what the Sussex PHP would have written: same checks, same format, same {id}_{task}.json
names. credit_report_gfactor.py reads them with --no-download --json-dir <out-dir>.

Usage:
  python3 import_datapipe.py ~/Downloads/gfactor --out-dir ~/gfactor_checked            # check into a folder
  python3 import_datapipe.py ~/Downloads/gfactor --out-dir ~/gfactor_checked --upload   # also copy new files to Sussex

The input is a folder of files downloaded from the DataPipe experiment's storage: loose .json files and/or the zip
archives DataPipe makes on Zenodo (it bundles older sessions once a record nears Zenodo's 100-file limit). Subfolders
are fine, and a file that appears both loose and in a zip is only checked once. Needs PHP on the PATH (e.g. from
Herd); --upload also needs SSH key access to unix.sussex.ac.uk.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from collections import defaultdict
from pathlib import Path

SSH_HOST = 'mel29@unix.sussex.ac.uk'
REMOTE_DIR = '/its/home/mel29/gfactor/tasks'
HERE = Path(__file__).resolve().parent
METHODS_DIR = HERE.parent / 'methods'
RUN_PHP = HERE / 'run_php.php'
# Same IDs the PHP scripts accept (5-6 digit SONA/random IDs or 24-character Prolific IDs)
FILE_RE = re.compile(r'^([0-9]{5,6}|[a-z0-9]{24})_(survey|span|dots|gabor|breath)\.json$')


def run_php(task, file_stem, exp_data, out_dir):
    """Run methods/<task>/<task>.php on the body the browser would have POSTed to it.

    Returns (path of the file the PHP wrote, None) or (None, reason it was rejected).
    """
    env = {'PATH': os.environ.get('PATH', '')}  # PHP sees these as $_SERVER - keep the shell's other variables out of it
    if task == 'survey' and isinstance(exp_data, dict):
        # save.js adds user_agent for DataPipe; survey.php rejects unknown fields and adds it from the request headers
        env['GFACTOR_USER_AGENT'] = str(exp_data.pop('user_agent', ''))
    out = out_dir / f'{file_stem}.json'
    out.unlink(missing_ok=True)  # never mistake a file from an earlier run for this one
    body = json.dumps({'file_name': file_stem, 'exp_data': exp_data})
    result = subprocess.run(['php', str(RUN_PHP), str(METHODS_DIR / task / f'{task}.php'), str(out_dir)],
                            input=body, capture_output=True, text=True, env=env, cwd=out_dir)
    if result.returncode == 2:
        sys.exit(result.stderr.strip())  # run_php.php couldn't adapt the PHP file - fix that before importing anything
    if result.returncode != 0:
        return None, f'PHP error: {(result.stderr or result.stdout).strip()[:200]}'
    return (out, None) if out.exists() else (None, f'rejected by {task}.php')


def list_server_files():
    result = subprocess.run(['ssh', '-o', 'BatchMode=yes', SSH_HOST, f'ls {REMOTE_DIR}'],
                            capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f'Could not list {REMOTE_DIR} on {SSH_HOST}: {result.stderr.strip()}')
    return set(result.stdout.split())


def collect(datapipe_dir):
    """{file name: [contents, ...]} for every .json file in the download, loose or inside DataPipe's zip archives."""
    found = defaultdict(list)
    for path in sorted(datapipe_dir.rglob('*')):
        if not path.is_file():
            continue
        if path.suffix == '.json':
            found[path.name].append(path.read_bytes())
        elif path.suffix == '.zip':
            with zipfile.ZipFile(path) as archive:
                for member in archive.infolist():
                    name = Path(member.filename).name  # only the name, whatever folders the archive holds
                    if not member.is_dir() and name.endswith('.json'):
                        found[name].append(archive.read(member))
    return found


def check_files(datapipe_dir, out_dir):
    """Run every study file in a DataPipe download through its task's PHP into out_dir.
    Returns (paths written, [(name, reason rejected)], [names ignored])."""
    checked, rejected, ignored = [], [], []
    for name, copies in sorted(collect(datapipe_dir).items()):
        match = FILE_RE.match(name)
        if not match:
            ignored.append(name)  # e.g. DataPipe metadata files
            continue
        if len(set(copies)) > 1:
            rejected.append((name, f'{len(copies)} different versions in the download - compare them by hand'))
            continue
        try:
            exp_data = json.loads(copies[0].decode('utf-8'))
        except (json.JSONDecodeError, UnicodeDecodeError) as e:
            rejected.append((name, f'not valid JSON: {e}'))
            continue
        path, reason = run_php(match[2], name[:-len('.json')], exp_data, out_dir)
        if reason:
            rejected.append((name, reason))
        else:
            checked.append(path)
    return checked, rejected, ignored


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('datapipe_dir', type=Path, help='folder of files downloaded from Zenodo or Google Drive')
    parser.add_argument('--out-dir', type=Path, help='put the checked files here (default: a temporary folder)')
    parser.add_argument('--upload', action='store_true', help='also copy checked files not yet on Sussex to the server')
    args = parser.parse_args()

    if not args.datapipe_dir.is_dir():
        sys.exit(f'Not a folder: {args.datapipe_dir}')
    if not shutil.which('php'):
        sys.exit("PHP isn't on the PATH - it's needed to run the study's own checks (methods/*/*.php)")
    out_dir = (args.out_dir or Path(tempfile.mkdtemp(prefix='gfactor_import_'))).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    checked, rejected, ignored = check_files(args.datapipe_dir, out_dir)

    print(f'{len(checked)} file(s) passed the checks and are in {out_dir}')
    if rejected:
        print(f'\n{len(rejected)} rejected:')
        print('\n'.join(f'  {name}: {reason}' for name, reason in rejected))
    if ignored:
        print(f'\nIgnored {len(ignored)} file(s) not named {{id}}_{{task}}.json: {", ".join(ignored)}')

    if not args.upload:
        print(f'\nFor credits: credit_report_gfactor.py --no-download --json-dir {out_dir}')
        return
    on_server = list_server_files()
    new = [path for path in checked if path.name not in on_server]
    already = sorted(path.name for path in checked if path.name in on_server)
    if already:
        print(f'\n{len(already)} already on the server, not overwritten (compare by hand if needed):')
        print('  ' + '\n  '.join(already))
    if new:
        subprocess.run(['scp', '-q', '-o', 'BatchMode=yes', *map(str, new), f'{SSH_HOST}:{REMOTE_DIR}/'], check=True)
    print(f'\nUploaded {len(new)} file(s) to {SSH_HOST}:{REMOTE_DIR}')


if __name__ == '__main__':
    main()
