import argparse
import os
import os.path
import shutil
import subprocess


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Wrapper to run rmats --task post'))
    parser.add_argument('--tmp-dir',
                        help='path to use for rmats --tmp', required=True)
    parser.add_argument('--out-dir',
                        help='path to use for rmats --od', required=True)
    parser.add_argument('--b1',
                        help='path to use for rmats --b1', required=True)
    parser.add_argument('--b2',
                        help='path to use for rmats --b2', required=True)
    parser.add_argument('--gtf',
                        help='path to use for rmats --gtf', required=True)
    parser.add_argument('--read-length',
                        help='value for rmats --readLength', type=int, required=True)
    parser.add_argument('--num-threads',
                        help='whether to use -t paired', type=int, required=True)
    parser.add_argument('--dot-rmats-files',
                        help='.rmats files to copy to --tmp-dir', required=True)
    parser.add_argument('--b1-bams',
                        help='bams to write to --b1', required=True)
    parser.add_argument('--b2-bams',
                        help='bams to write to --b2', required=False)

    args = parser.parse_args()
    args.dot_rmats_files = args.dot_rmats_files.split(',')
    args.b1_bams = args.b1_bams.split(',')
    if args.b2_bams:
        args.b2_bams = args.b2_bams.split(',')

    return args


def write_b1_or_b2(bams, file_name):
    abs_paths = list()
    for bam in bams:
        abs_path = os.path.abspath(bam)
        abs_paths.append(abs_path)

    comma_paths = ','.join(abs_paths)
    with open(file_name, 'wt') as f_handle:
        f_handle.write('{}\n'.format(comma_paths))


def run_rmats_post(args):
    write_b1_or_b2(args.b1_bams, args.b1)
    if args.b2_bams:
        write_b1_or_b2(args.b2_bams, args.b2)
    else:
        with open(args.b2, 'wt'):
            pass  # create empty file

    copy_dot_rmats(args.tmp_dir, args.dot_rmats_files)
    rmats_command = ['rmats.py', '--tmp', args.tmp_dir, '--od', args.out_dir,
                     '--b1', args.b1, '--task', 'post', '--readLength', str(args.read_length),
                     '--gtf', args.gtf, '--variable-read-length', '--novelSS', '--nthread',
                     str(args.num_threads)]
    if args.b2_bams:
        rmats_command.extend(['--b2', args.b2])
    else:
        rmats_command.append('--statoff')

    print(rmats_command)
    subprocess.run(rmats_command, check=True)


def copy_dot_rmats(tmp_dir, dot_rmats_files):
    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)

    for path in dot_rmats_files:
        print('cp {} {}'.format(path, tmp_dir))
        shutil.copy(path, tmp_dir)


def main():
    args = parse_args()
    run_rmats_post(args)


if __name__ == '__main__':
    main()
