"""
Check g-factor study files collected through DataPipe (Google Drive) and upload them to the Sussex server.

While the study runs from GitHub Pages, each task saves to DataPipe instead of the Sussex PHP scripts (see
methods/save.js). This script runs every file through that task's original PHP script via run_php.php, so what ends
up in /its/home/mel29/gfactor/tasks is exactly what the PHP would have written: same checks, same format, same
{id}_{task}.json names. credit_report_gfactor.py and the analysis need no changes.

Usage:
  python3 import_datapipe.py ~/Downloads/gfactor             # check only: reports what would be uploaded
  python3 import_datapipe.py ~/Downloads/gfactor --upload    # also copies the new files to the server

The input is the DataPipe experiment folder downloaded from Google Drive and unzipped (subfolders are fine).
Needs PHP on the PATH (e.g. from Herd) and SSH key access to unix.sussex.ac.uk.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
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


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('datapipe_dir', type=Path, help='unzipped DataPipe folder downloaded from Google Drive')
    parser.add_argument('--upload', action='store_true', help='copy the checked, new files to the Sussex server')
    parser.add_argument('--out-dir', type=Path, help='keep the checked files here (default: a temporary folder)')
    args = parser.parse_args()

    if not args.datapipe_dir.is_dir():
        sys.exit(f'Not a folder: {args.datapipe_dir}')
    if not shutil.which('php'):
        sys.exit("PHP isn't on the PATH - it's needed to run the study's own checks (methods/*/*.php)")
    out_dir = (args.out_dir or Path(tempfile.mkdtemp(prefix='gfactor_import_'))).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    on_server = list_server_files()

    new, already_on_server, rejected, ignored, seen = [], [], [], [], set()
    for path in sorted(args.datapipe_dir.rglob('*.json')):
        match = FILE_RE.match(path.name)
        if not match:
            ignored.append(path.name)  # e.g. DataPipe metadata files
            continue
        if path.name in seen:
            rejected.append((path.name, f'duplicate name ({path})'))
            continue
        seen.add(path.name)
        try:
            exp_data = json.loads(path.read_text(encoding='utf-8'))
        except (json.JSONDecodeError, UnicodeDecodeError) as e:
            rejected.append((path.name, f'not valid JSON: {e}'))
            continue
        checked, reason = run_php(match[2], path.stem, exp_data, out_dir)
        if reason:
            rejected.append((path.name, reason))
        elif path.name in on_server:
            already_on_server.append(path.name)
        else:
            new.append(checked)

    print(f'Checked files are in {out_dir}')
    print(f'{len(new)} new file(s) passed the checks')
    if already_on_server:
        print(f'\n{len(already_on_server)} already on the server, not overwritten (compare by hand if needed):')
        print('  ' + '\n  '.join(already_on_server))
    if rejected:
        print(f'\n{len(rejected)} rejected (not uploaded):')
        print('\n'.join(f'  {name}: {reason}' for name, reason in rejected))
    if ignored:
        print(f'\nIgnored {len(ignored)} file(s) not named {{id}}_{{task}}.json: {", ".join(ignored)}')

    if not args.upload:
        if new:
            print('\nCheck only - run again with --upload to copy the new files to the server.')
        return
    if new:
        subprocess.run(['scp', '-q', '-o', 'BatchMode=yes', *map(str, new), f'{SSH_HOST}:{REMOTE_DIR}/'], check=True)
        print(f'\nUploaded {len(new)} file(s) to {SSH_HOST}:{REMOTE_DIR}')


if __name__ == '__main__':
    main()
