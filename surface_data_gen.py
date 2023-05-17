import sys
from utils.masif.generate_prot_ply import compute_inp_surface
import os.path as osp
from glob import glob
import argparse


parser = argparse.ArgumentParser()

parser.add_argument(
    '--pkt_file', action='store',required=False,type=str,default='./example/1bxm/1bxm_protein.pdb',
    help='protein pdb file'
)
parser.add_argument(
    '--sdf_file', action='store',required=False,type=str,default='./example/1bxm/1bxm_protein.pdb',
    help='ligand sdf file'
)
parser.add_argument(
    '--out_dir', action='store',required=False,type=str,default='./',
    help='output directory where the surface file will be generated'
)
parser.add_argument(
    '--out_name', action='store',required=False,type=str,default='surface',
    help='out surface name ./outdir/{surface}.ply'
)


args = parser.parse_args()
compute_inp_surface(args.pkt_file, args.sdf_file, args.out_dir, \
    out_name=args.out_name, mesh_res=1.5,dist_threshold=8.0)