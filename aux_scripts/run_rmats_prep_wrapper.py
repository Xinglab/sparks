import argparse
import os
import os.path
import shutil
import subprocess


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Wrapper to run rmats --task prep'))
    parser.add_argument(
        '--bam',
        help='path to .bam input file', required=True)
    parser.add_argument('--dot-rmats',
                        help='path for output .rmats file', required=True)
    parser.add_argument('--tmp-dir',
                        help='path to use for rmats --tmp', required=True)
    parser.add_argument('--out-dir',
                        help='path to use for rmats --od', required=True)
    parser.add_argument('--b1',
                        help='path to use for rmats --b1', required=True)
    parser.add_argument('--gtf',
                        help='path to use for rmats --gtf', required=True)
    parser.add_argument('--read-length',
                        help='value for rmats --readLength', type=int, required=True)
    parser.add_argument('--is-paired',
                        help='whether to use -t paired', action='store_true', required=True)

    return parser.parse_args()


def write_b1(bam, file_name):
    abs_path = os.path.abspath(bam)
    with open(file_name, 'wt') as f_handle:
        f_handle.write('{}\n'.format(abs_path))


def run_rmats_prep(args):
    write_b1(args.bam, args.b1)
    rmats_command = ['rmats.py', '--tmp', args.tmp_dir, '--od', args.out_dir,
                     '--b1', args.b1, '--task', 'prep', '--readLength', str(args.read_length),
                     '--gtf', args.gtf, '--variable-read-length', '--novelSS']
    if args.is_paired:
        rmats_command.extend(['-t', 'paired'])
    else:
        rmats_command.extend(['-t', 'single'])

    print(rmats_command)
    subprocess.run(rmats_command, check=True)

    print('rm {}'.format(args.b1))
    os.remove(args.b1)


def move_dot_rmats(tmp_dir, final_dot_rmats_path):
    dot_rmats_path = None
    tmp_output_files = os.listdir(tmp_dir)
    dot_rmats_names = [x for x in tmp_output_files if x.endswith('.rmats')]
    if len(dot_rmats_names) != 1:
        raise Exception('expected 1 .rmats output file: {}'.format(dot_rmats_names))

    dot_rmats_path = os.path.join(tmp_dir, dot_rmats_names[0])
    print('mv {} {}'.format(dot_rmats_path, final_dot_rmats_path))
    shutil.move(dot_rmats_path, final_dot_rmats_path)


def main():
    args = parse_args()
    run_rmats_prep(args)
    move_dot_rmats(args.tmp_dir, args.dot_rmats)


if __name__ == '__main__':
    main()
