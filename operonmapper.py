"""
Use OperonMapper from your CLI
https://biocomputo.ibt.unam.mx/operon_mapper/
"""
import os
import io
import re
import httpx
import time
import tarfile


OPERON_MAPPER_URL = 'https://biocomputo.ibt.unam.mx/operon_mapper'


def _submit(
    fastafile: io.TextIOWrapper,
    gfffile: io.TextIOWrapper = None,
    description: str = 'An Operon Search',
    email: str = None

):
    form_action = f'{OPERON_MAPPER_URL}/capta_forma_01.pl'

    data = {
        'fastaseq': fastafile.read(),
        'descri': description,
        'todoscomprimido': 'si'
    }

    if gfffile:
        data['gffseq'] = gfffile.read()
    if email:
        data['email1'] = email

    headers = {
        'Origin': 'https://biocomputo.ibt.unam.mx',
        'Referer': 'https://biocomputo.ibt.unam.mx/operon_mapper/',
    }

    r = httpx.post(url=form_action, headers=headers, data=data)
    r.raise_for_status()

    if match := re.search(r'out_(?P<id>\d+)\.html', r.text, re.MULTILINE):
        return match.groupdict()['id']

    raise Exception('Link to results page not found\n' + r.text)


def _read(out_id, out_dir):
    """
    Out id example http://biocomputo.ibt.unam.mx/operon_mapper/out/2614089.tar.gz
    """
    os.makedirs(out_dir, exist_ok=True)
    rsp = httpx.get(f'{OPERON_MAPPER_URL}/out/{out_id}.tar.gz')
    rsp.raise_for_status()

    with tarfile.open(fileobj=io.BytesIO(rsp.read())) as targz:
        targz.extractall(out_dir)

    out_folder = f'{out_dir}/{out_id}'
    files = os.listdir(out_folder)

    for file in files:
        fullname = out_folder + '/' + file
        if os.path.isdir(fullname):
            continue
        os.rename(fullname, out_dir + '/' + file.rsplit('_', maxsplit=1)[0])

    os.rmdir(out_dir + '/' + out_id)

    return sorted(file for file in os.listdir(out_dir))


def _peek(out_id):
    """ 
    Checks wether results is ready or not 
    recommended peek time - 1/30sec
    """
    rsp = httpx.get(f'{OPERON_MAPPER_URL}/out/out_{out_id}.html')
    rsp.raise_for_status()
    text = rsp.text

    if 'Output files' in text and 'Compressed file with all the above' in text:
        return True

    return False


def _find_operons(args):
    if not args.description:
        args.description = os.path.basename(args.fasta)

    fastafile = args.fasta and open(args.fasta)
    gfffile = args.gff and open(args.gff)

    args.task_id = _submit(fastafile, gfffile, args.description, args.email)

    fastafile and fastafile.close()
    gfffile and gfffile.close()

    print('Operon finder task submitted')
    print('\tTask id:', args.task_id)
    print('Visit', f'{OPERON_MAPPER_URL}/out/out_{args.task_id}.html', 'to manually check status')

    return _do_pooling(args)


def _prettysleep(n, freq=3):
    rate = 1/freq
    for _ in range(n*freq):
        time.sleep(rate)
        yield rate


def _do_pooling(args):

    if not args.output:
        args.output = 'operonmapper_output_' + args.task_id

    task_id = args.task_id
    t = 0
    while not _peek(task_id):
        for dt in _prettysleep(15):
            t += dt
            print(f'\rResults are not ready yet, keep waiting ({t:.2f}sec passed)', end='')

    print()

    files = _read(task_id, args.output)

    print(f'Result is ready, it took only {t} seconds\n')
    print('Output files are:')
    for f in files:
        print('\t', args.output + '/' + f)


def main():
    import argparse

    main_parser = argparse.ArgumentParser()

    subs = main_parser.add_subparsers(required=True)
    parser = subs.add_parser('start')
    parser.add_argument('fasta', help='path to fasta file')
    parser.add_argument('--gff', help='path to gff file')
    parser.add_argument('-d', '--description',
                        help='opperon mapper job name (defaults to fasta file name)')
    parser.add_argument('--email', default=None,
                        help='email to be notified regarding job status')
    parser.add_argument('-o', '--output', default=None,
                        help='output directory (defaults to task_id)')


    parser.set_defaults(func=_find_operons)

    continue_parser = subs.add_parser('continue')
    continue_parser.add_argument('task_id', help='task id to continue pooling')
    continue_parser.set_defaults(func=_do_pooling)
    continue_parser.add_argument('-o', '--output', default=None,
                                 help='(defaults to task_id)')

    args = main_parser.parse_args()

    return args.func(args)


if __name__ == '__main__':
    main()
