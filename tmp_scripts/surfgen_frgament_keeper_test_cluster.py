# on cluster
# python -u surfgen_frgament_keeper_test.py --check_point ./ckpt_frag/crossdock_val_159.pt 
import os
import argparse
from glob import glob
from easydict import EasyDict
from Bio import BiopythonWarning
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Selection import unfold_entities
from rdkit import Chem
import torch
# from feats.protein import get_protein_feature_v2
from Bio.PDB import NeighborSearch, Selection
from utils.protein_ligand import parse_rdmol, parse_sdf_file
from utils.data import torchify_dict, ProteinLigandData
from copy import deepcopy
import shutil
import numpy as np
from tqdm.auto import tqdm
from utils.transforms import *
from utils.misc import load_config
from utils.reconstruct import *
from utils.datasets.surfdata import SurfGenDataset
from models.surfgen import SurfGen

from utils.sample import get_init, get_next, logp_to_rank_prob
from utils.sample import STATUS_FINISHED, STATUS_RUNNING
import os.path as osp
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
import warnings

import pickle
def write_pkl(list,file):
    with open(file,'wb') as f:
        pickle.dump(list,f)
        print('pkl file saved at {}'.format(file))
def read_pkl(file):
    with open(file,'rb') as f:
        data = pickle.load(f)
    return data

from plyfile import PlyData
PlyData

def read_sdf(sdf_file):
    supp = Chem.SDMolSupplier(sdf_file)
    mols_list = [i for i in supp]
    return mols_list

def write_sdf(mol_list,file):
    writer = Chem.SDWriter(file)
    for i in mol_list:
        writer.write(i)
    writer.close()

def pdb_to_pocket_data(ply_file, keep_frag):
    '''
    use the sdf_file as the center 
    '''
    protein_dict = read_ply(ply_file)
    keep_frag_mol = read_sdf(keep_frag)[0]
    ligand_dict = parse_rdmol(keep_frag_mol)
    data = ProteinLigandData.from_protein_ligand_dicts(
        protein_dict = torchify_dict(protein_dict),
        ligand_dict = torchify_dict(ligand_dict)
    )
    return data

parser = argparse.ArgumentParser()

parser.add_argument('--config', type=str, default='./configs/sample.yml')
parser.add_argument('--outdir', type=str, default='./outputs')
parser.add_argument('--device', type=str, default='cuda')
parser.add_argument('--check_point',type=str,default='./ckpt_frag/val_159.pt') #22, 53, 119, 235,linker_val_37.pt
parser.add_argument('--SDF_dirname', type=str,default='SDF')
parser.add_argument('--sdf_filename', type=str,default='gen')
parser.add_argument('--surf_path', type=str,default='./data/bindingmoad_surface_8')
parser.add_argument('--frag_path', type=str,default='./data/moad_linkers/linker_fragments')
parser.add_argument('--index_data', type=str,default='./data/moad_linkers/index_data.pkl')
parser.add_argument('--index_file', type=str,default='./data/moad_linkers/linker_index.pkl')
parser.add_argument('--frag_freeze', action='store_true')

args = parser.parse_args()

# load configs, utils 
config = load_config(args.config)
ckpt = torch.load(args.check_point, map_location=args.device)
config_train = ckpt['config']
protein_featurizer = FeaturizeProteinAtom()
ligand_featurizer = FeaturizeLigandAtom()
composer = AtomComposer(protein_featurizer.feature_dim, ligand_featurizer.feature_dim, config_train.model.encoder.knn)
transform = Compose([
    RefineData(),
    LigandCountNeighbors(),
    ligand_featurizer,
    protein_featurizer
])
mask = LigandMaskZero()
masking = Compose([
    mask, 
    composer
])
# model loading
model = SurfGen(
    ckpt['config'].model, 
    num_classes = 7,
    num_bond_types = 3,
    protein_atom_feature_dim = protein_featurizer.feature_dim,
    ligand_atom_feature_dim = ligand_featurizer.feature_dim,
).to(args.device)
model.load_state_dict(ckpt['model'])
print('Num of parameters is {0:.4}M'.format(np.sum([p.numel() for p in model.parameters()]) /100000 ))


frag_freeze = args.frag_freeze 
# . the data
index_data = read_pkl(args.index_data)
for data in index_data:
    surf_fn = data.surf_fn
    frag_fn = data.frag_fn
    lig_fn = data.lig_fn
    surf_file = osp.join(args.surf_path, surf_fn)
    frag_file = osp.join(args.frag_path, frag_fn)
    
    # you can iterate the data for conditional generation on different targets, but we use the default 0
    if args.frag_freeze:
        print('Fragment Freezing!')
        # freeze the part of fragment
        mol = read_sdf(frag_file)[0]
        all_id = np.arange(mol.GetNumAtoms())
        anchor_idx = [float(i) for i in mol.GetProp('anchor_idx').split('_')]
        for anchor_id in anchor_idx:
            all_id = all_id[~(all_id==anchor_id)]
        freeze = torch.tensor(all_id).to(args.device)
    else:
        # no freeze for fragment
        freeze = None
    try:
        data = transform(data)
        data = transform_data(data, masking)
    except:
        continue
        print('transform data failed')
    # generation
    np.seterr(invalid='ignore') 
    pool = EasyDict({
        'queue': [],
        'failed': [],
        'finished': [],
        'duplicate': [],
        'smiles': set(),
    })

    print('Start to generate!')
    init_data_list = get_next(
                data.to(args.device), 
                model = model,
                transform = composer,
                threshold = config.sample.threshold,
                frontier_threshold=0.0,
                freeze = freeze
            )
    pool.queue = init_data_list
    #rint('Start to generate novel molecules with 3D conformation located in the protein pocket!')
    #print('The protein pocket is {}, init length is {}'.format(data.protein_filename, len(init_data_list)))
    global_step = 0 
    while len(pool.finished) < config.sample.num_samples:
        global_step += 1
        if global_step > config.sample.max_steps:
            break
        queue_size = len(pool.queue)
        # # sample candidate new mols from each parent mol
        queue_tmp = []
        for data in pool.queue:
            nexts = []
            data_next_list = get_next(
                data.to(args.device), 
                model = model,
                transform = composer,
                threshold = config.sample.threshold,
                freeze = freeze
            )

            for data_next in data_next_list:
                if data_next.status == STATUS_FINISHED:
                    try:
                        rdmol = reconstruct_from_generated_with_edges(data_next)
                        data_next.rdmol = rdmol
                        mol = Chem.MolFromSmiles(Chem.MolToSmiles(rdmol))
                        smiles = Chem.MolToSmiles(mol)
                        data_next.smiles = smiles
                        if smiles in pool.smiles:
                            #print('Duplicate molecule: %s' % smiles)
                            pool.duplicate.append(data_next)
                        elif '.' in smiles:
                            print('Failed molecule: %s' % smiles)
                            pool.failed.append(data_next)
                        else:   # Pass checks
                            print('Success: %s' % smiles)
                            pool.finished.append(data_next)
                            pool.smiles.add(smiles)
                    except MolReconsError:
                        #print('Reconstruction error encountered.')
                        pool.failed.append(data_next)
                elif data_next.status == STATUS_RUNNING:
                    nexts.append(data_next)

            queue_tmp += nexts
        prob = logp_to_rank_prob(np.array([p.average_logp[2:] for p in queue_tmp]),)  # (logp_focal, logpdf_pos), logp_element, logp_hasatom, logp_bond
        n_tmp = len(queue_tmp)
        if n_tmp == 0:
            print('failures!')
            break
        else:
            next_idx = np.random.choice(np.arange(n_tmp), p=prob, size=min(config.sample.beam_size, n_tmp), replace=False)
        pool.queue = [queue_tmp[idx] for idx in next_idx]
    try:
        ckpt_name = args.check_point.split('/')[-1][:-3]
        if args.frag_freeze:
            ckpt_name = ckpt_name + '_kept_frag'
        out_dir = osp.join(args.outdir,frag_fn)
        os.makedirs(out_dir, exist_ok=True)
        sdf_name = frag_fn[:-3] + f'_{ckpt_name}_.sdf'
        sdf_file = os.path.join(out_dir,sdf_name)
        writer = Chem.SDWriter(sdf_file)
        for j in range(len(pool['finished'])):
            writer.write(pool['finished'][j].rdmol)
        writer.close()

        SDF_dir = os.path.join(out_dir, ckpt_name) 
        os.makedirs(SDF_dir, exist_ok=True)
        for j in range(len(pool['finished'])):
            writer = Chem.SDWriter(os.path.join(SDF_dir,f'{j}.sdf'))
            writer.write(pool['finished'][j].rdmol)
        writer.close()
    except:
        print('write the generated mols failed')
