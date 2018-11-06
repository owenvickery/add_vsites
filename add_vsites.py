import os, sys 
import numpy as np
import argparse
from subprocess import Popen, PIPE
import subprocess, shlex
from time import gmtime, strftime
from shutil import copyfile
from pathlib import Path
from collections import Counter  
import re
import time

convert=dict([('POPC', ['NC3', 'N  ']), ('POPG', ['CH3', 'C13'])])
lipid_residues=dict([('POPE', 131), ('POPG',131), ('POPC',144), ('CHL',84), ('PIP2',148), ('CARD', 248), ('UDP2',178), ('DMPC', 128), ('DPPC', 140)])
solvent=dict([('SOL',3),('NA',1), ('k',1), ('CL',1)])
cations=dict([('NA',1), ('k',1)])
anions=dict([('CL',1)])
amino_acids = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'NME': 'NME', 'ACE': 'ACE' }

# def gromacs(cmd):
# 	print('\nrunning gromacs: \n '+cmd)
# 	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# 	err, out = output.communicate()
# 	exitcode = output.returncode
# 	out=out.decode("utf-8")
# 	checks = open('gromacs_outputs', 'a')
# 	checks.write(out)

def make_min():
	if not os.path.exists('em.mdp'):
		with open('em.mdp', 'w') as em:
			em.write('integrator = steep\nnsteps     = 5000\nemtol      = 1000\nemstep     = 0.001')

def read_topology(itp):
	if Path(itp+'.top').exists():
		read=False
		with open(itp+'.itp', 'w') as itp_write:
			for line in open(itp+'.top', 'r').readlines():
				if len(line.split()) > 1: 
					if read == False and line.split()[1] == 'moleculetype':
						read = True
					if read == True and line.split()[1] == 'system':
						read = False
					if line.split()[0] == 'Other':
							line= re.sub('Other', 'ligand',line)
				if read ==True:
					itp_write.write(line)

def read_pdb(pdb_file, run):
	box, ace, nme=False, False,False
	pdb=[[],[],[],[]]
	lipid_order, solvent_order=[],[]
	prot,lig,lip, sol={},{},{},{}
	box_line=''
	for line in open(pdb_file, 'r').readlines():
		if not line[0] in ['#', '@']:
			if line.split()[0] == 'CRYST1' and box == False:
				box_line=line
				box=True
			if line.split()[0] == 'ATOM':
				if line.split()[3][:4] in amino_acids:
					if line.split()[3][:4] =='ACE':
						ace=True
					if line.split()[3][:4] =='NME':
						nme=True
					if run == 1:
						if 'OT' in line:
							line= re.sub(r'OT[0-9]', 'O  ',line)				
					pdb[0].append(line)
				elif line.split()[3][:4] == args.ligand:
					pdb[1].append(line)
				elif line.split()[3][:4] in lipid_residues:
					if line.split()[3][:4] in convert and convert[line.split()[3][:4]][0] in line:
						line= re.sub(convert[line.split()[3][:4]][0], convert[line.split()[3][:4]][1],line)
					pdb[2].append(line)
					if run == 2:
						if not line.split()[3][:4] in lip:
							lip[line.split()[3][:4]] = 1
							lipid_order.append(line.split()[3][:4])
						else:
							lip[line.split()[3][:4]] += 1
				elif line.split()[3][:4] in solvent:
					pdb[3].append(line)
					if run == 2:
						if not line.split()[3][:4] in sol:
							sol[line.split()[3][:4]] = 1
							solvent_order.append(line.split()[3][:4])
						else:
							sol[line.split()[3][:4]] += 1
	return box_line,pdb, ace, nme,prot,lig, lip, sol, lipid_order,solvent_order

def pdb2gmx():
	box_line,pdb, ace, nme,prot,lig, lip, sol, lipid_order,solvent_order = read_pdb(args.f, 1)
	structure=['protein','ligand','lipids','solvent']
	for val, section in enumerate(pdb):
		if len(section)>0:
			with open(structure[val]+'.pdb', 'w') as file_write:
				file_write.write('TITLE     '+structure[val]+'\n'+box_line)
				for j in section:
					file_write.write(j)
			
		if structure[val] in ['ligand','lipids']:
			os.system('gmx pdb2gmx -vsite h -f '+structure[val]+'.pdb -o '+structure[val]+'_vs.pdb -water none -p '+structure[val]+'_topology.top << EOF \n1\nEOF')
		cter, nter='1','1'
		if ace ==True:
			nter='3'
		if args.nter==True:
			nter = '0'
		if nme ==True:
			cter='4'
		if args.cter==True:
			cter = '0'		 	
		if structure[val]== 'protein':
			os.system('pwd')
			os.system('gmx pdb2gmx -ter -ignh -vsite h -f '+structure[val]+'.pdb -o '+structure[val]+'_vs.pdb -water none -p '+structure[val]+'_topology.top << EOF \n1\n'+nter+'\n'+cter+'\nEOF')
	os.system('rm *itp')


def make_system():
	os.system('cat protein_vs.pdb ligand_vs.pdb lipids_vs.pdb solvent.pdb > system_vs_concat.pdb')
	box_line,pdb, ace, nme,prot,lig, lip, sol, lipid_order, solvent_order  = read_pdb('system_vs_concat.pdb', 2)
	with open('system_vs_final.pdb', 'w') as system_write:
		system_write.write('TITLE     system with virtual sites\n'+box_line)
		for val, section in enumerate(pdb):
			if len(section)>0:
				for j in section:
					system_write.write(j)
	for top in ['protein_topology','ligand_topology']:
		read_topology(top)
	with open('topol.top', 'w') as topol_write:
		topol_write.write('#include \"charmm36-jul2017-updated.ff/forcefield.itp\"\n#include \"charmm36-jul2017-updated.ff/lipids_vs.itp\"\n')
		if Path('protein_topology.itp').exists():
			topol_write.write('#include \"protein_topology.itp\"\n')
		if Path('ligand_topology.itp').exists():
			topol_write.write('#include \"ligand_topology.itp\"\n')	
		topol_write.write('#include \"charmm36-jul2017-updated.ff/tip3p.itp\"\n\n#ifdef POSRES_WATER\n; Position restraint for each water oxygen\n[ position_restraints ]\n\
;  i funct       fcx        fcy        fcz\n   1    1       1000       1000       1000\n#endif\n\n#include \"charmm36-jul2017-updated.ff/ions.itp\"\n\n\
[ system ]\n; Name\nGreat Red Oystrich Makes All Chemists Sane\n\n[ molecules ]\n; Compound        #mols\n')
		lip_number = {k: lip[k]/lipid_residues[k] for k in lip.keys() & lipid_residues}
		sol_number = {k: sol[k]/solvent[k] for k in sol.keys() & solvent}
		if len(pdb[0]) > 0:
			topol_write.write('protein\t1\n')
		if len(pdb[1]) > 0:
			topol_write.write('ligand\t1\n')	
		for i in lipid_order:
			topol_write.write(i+'\t'+str(int(lip_number[i]))+'\n')
		for i in solvent_order:
			topol_write.write(i+'\t'+str(int(sol_number[i]))+'\n')	
	minimise()

def minimise():
	make_min()
	try: 
		os.makedirs('minimisation')
	except:
		print('minimisation')
	os.system('gmx grompp -p topol.top -c system_vs_final -o minimisation/system_vs_minimised -f em.mdp')
	os.chdir('minimisation')
	os.system('gmx mdrun -v -deffnm system_vs_minimised')


###### forcefield

def input(files, section, idih):
	count=False
	d=[]
	moltype, count=False, False
	inp=open(files)
	for line in open(files, 'r').readlines():
		if len(line.split()) >=2:
				if line.split()[1] == section:
						moltype=True
				if moltype==True and line.split()[0][0]!='[' and idih==False:
					if line.split()[0][0] not in (';', '#', '[',''):
						d.append(line)
		if moltype==True and len(line.split())==0 and idih==False:
			return d
		if moltype==True and len(line.split())==0 and idih==True:
			moltype=True
			idih=False
	return d

def forceread(files):
	d=[]
	for line in open(files, 'r').readlines():
		d.append(line)
	return d

def merged_add(forcefield, to_add):
	print('adding to merged.rtp')
	merged_rtp=forceread(forcefield+to_add)
	new_merged = open(to_add, 'w')
	already_exists=False
	for line_ori_merged in merged_rtp:
		if len(line_ori_merged.split())==3:
				if moleculetype[0].split()[0] == line_ori_merged.split()[1]:
					already_exists=True
		new_merged.write(line_ori_merged)
	if already_exists == False:
		new_merged.write('\n\n[ '+moleculetype[0].split()[0]+' ]\n')
		new_merged.write('  [ atoms ]\n')
		for atom_line in atoms:
			line=atom_line.split()[4]+'\t\t'+atom_line.split()[1]+'\t'+atom_line.split()[6]+'\t'+atom_line.split()[0]+'\t'+atom_line.split()[8]+'\t'+atom_line.split()[9]+'\t'+atom_line.split()[10]+'\n'
			new_merged.write(line)
		new_merged.write('  [ bonds ]\n')
		for bond_line in bonds:
			line=atoms[int(bond_line.split()[0])-1].split()[4]+'\t'+atoms[int(bond_line.split()[1])-1].split()[4]+'\n'
			new_merged.write(line)
		if len(idihtypes_lig)>0:
			new_merged.write('  [ impropers ]\n')
			for idih_line in idihtypes_lig:
				line=atoms[int(idih_line.split()[0])-1].split()[4]+'\t'+atoms[int(idih_line.split()[1])-1].split()[4]+'\t'+atoms[int(idih_line.split()[2])-1].split()[4]+'\t'+atoms[int(idih_line.split()[3])-1].split()[4]+'\n'
				new_merged.write(line)

def nonbonded_add(forcefield, to_add):
	print('adding to ffnonbonded.itp')
	atomtypes=input(charmm, 'atomtypes', False)
	pairtypes=input(charmm, 'pairtypes', False)
	nonbonded_itp=forceread(forcefield+to_add)
	new_nonbonded = open(to_add, 'w')
	atom_types, pair_types, cont=False,False, False
	check=[]
	old_atomtypes,old_pairtypes=[],[]
	for nonbonded_line in nonbonded_itp:
		if nonbonded_line.replace('\n', '')=='[ atomtypes ]':
			atom_types=True
		elif atom_types ==False and pair_types==False:
			new_nonbonded.write(nonbonded_line)
		elif nonbonded_line.replace('\n', '')=='[ pairtypes ]':
			new_nonbonded.write('\n')
			atom_types, pair_types=False,True
		if atom_types==True and len(nonbonded_line.replace('\n', ''))!=0:
			new_nonbonded.write(nonbonded_line)
			old_atomtypes.append(nonbonded_line.replace('\n', ''))

		if pair_types==True and len(nonbonded_line.replace('\n', ''))!=0:
			new_nonbonded.write(nonbonded_line)
			old_pairtypes.append(nonbonded_line.replace('\n', ''))

		if atom_types==True and len(nonbonded_line.replace('\n', ''))==0:
			for atomtype_line in atomtypes:
				for old_atomtypes_line in old_atomtypes:
					if atomtype_line.split()[0] == old_atomtypes_line.split()[0]:
						if "{:.3f}".format(float(atomtype_line.split()[5])) == "{:.3f}".format(float(old_atomtypes_line.split()[5])) and "{:.3f}".format(float(atomtype_line.split()[6])) == "{:.3f}".format(float(old_atomtypes_line.split()[6])):
							cont=True
						else:
							print('not copied:   \n',atomtype_line.replace('\n', ''), '\n', old_atomtypes_line)
							cont=True
				if cont==False:
					new_nonbonded.write(atomtype_line)
				cont=False
		if pair_types==True and len(nonbonded_line.replace('\n', ''))==0:
			for pairtype_line in pairtypes:
				if any(pairtype_line.split()[0] in s for s in old_pairtypes) or any(pairtype_line.split()[1] in s for s in old_pairtypes):

					for old_pairtypes_line in old_pairtypes:
						if pairtype_line.split()[0] == old_pairtypes_line.split()[0] and pairtype_line.split()[1] == old_pairtypes_line.split()[1]:
							if "{:.3f}".format(float(pairtype_line.split()[3])) == "{:.3f}".format(float(old_pairtypes_line.split()[3])) and "{:.3f}".format(float(pairtype_line.split()[4])) == "{:.3f}".format(float(old_pairtypes_line.split()[4])):
								cont=True
							else:
								print('not copied:   \n',pairtype_line.replace('\n', ''), '\n', old_pairtypes_line)
								cont=True
						if pairtype_line.split()[0] == old_pairtypes_line.split()[1] and pairtype_line.split()[1] == old_pairtypes_line.split()[0]:
							if "{:.3f}".format(float(pairtype_line.split()[3])) == "{:.3f}".format(float(old_pairtypes_line.split()[3])) and "{:.3f}".format(float(pairtype_line.split()[4])) == "{:.3f}".format(float(old_pairtypes_line.split()[4])):
								cont=True
							else:
								print('not copied:   \n',pairtype_line.replace('\n', ''), '\n', old_pairtypes_line)
								cont=True
				if cont==False:
					new_nonbonded.write(pairtype_line)
				cont=False
			break

def bonded_add(forcefield, to_add):
	print('adding to ffbonded.itp')
	bonded_itp=forceread(forcefield+to_add)
	new_bonded = open(to_add, 'w')
	bond_types,constraint_types, angle_types,dih_types, idih_types, idih, cont =False,False,False,False, False,False, False
	constr, angle, dih, idih,dih_check, test=False,False,False,False, False, False
	check_bond, check_angle, check_dih, check_idih=[],[],[],[]
	old_bondtypes, old_constraint, old_angle, old_dih, old_idih=[],[],[], [], []
	idih =0
	for bonded_line in bonded_itp:

		if bonded_line.replace('\n', '')=='[ bondtypes ]':
			new_bonded.write('\n')
			bond_types=True
		elif bond_types ==False and constraint_types==False and angle_types ==False and dih_types==False:
			new_bonded.write(bonded_line)
		elif bonded_line.replace('\n', '')=='[ constrainttypes ]':
			new_bonded.write('\n')
			bond_types, constraint_types=False,True
		elif bonded_line.replace('\n', '')=='[ angletypes ]':
			new_bonded.write('\n')
			constraint_types, angle_types=False,True
		elif bonded_line.replace('\n', '')=='[ dihedraltypes ]': #[ dihedraltypes ]dihedraltypes
			new_bonded.write('\n')
			angle_types, dih_types=False,True
			idih +=1
		if bond_types==True and len(bonded_line.replace('\n', ''))!=0:
			new_bonded.write(bonded_line)
			old_bondtypes.append(bonded_line.replace('\n', ''))
		if constraint_types==True and len(bonded_line.replace('\n', ''))!=0:
			new_bonded.write(bonded_line)
			old_constraint.append(bonded_line.replace('\n', ''))
		if angle_types==True and len(bonded_line.replace('\n', ''))!=0:
			new_bonded.write(bonded_line)
			old_angle.append(bonded_line.replace('\n', ''))
		if dih_types==True and len(bonded_line.replace('\n', ''))!=0 and idih==1:
			new_bonded.write(bonded_line)
			old_dih.append(bonded_line.replace('\n', ''))
		if dih_types==True and len(bonded_line.replace('\n', ''))!=0 and idih==2:
			new_bonded.write(bonded_line)
			old_idih.append(bonded_line.replace('\n', ''))



		if bond_types==True and len(bonded_line.replace('\n', ''))==0:
			for bondtype_line in bondtypes:
				if any(bondtype_line.split()[0] in s for s in old_bondtypes) or any(bondtype_line.split()[1] in s for s in old_bondtypes):
					for old_bondtypes_line in old_bondtypes:
						if old_bondtypes_line[0][0] not in (';', '#', '[',''):
							if bondtype_line.split()[0] == old_bondtypes_line.split()[0] and bondtype_line.split()[1] == old_bondtypes_line.split()[1]:
								if "{:.3f}".format(float(bondtype_line.split()[3])) == "{:.3f}".format(float(old_bondtypes_line.split()[3])) and "{:.3f}".format(float(bondtype_line.split()[4])) == "{:.3f}".format(float(old_bondtypes_line.split()[4])):
									cont=True
								else:
									print('not copied:   \n',bondtype_line.replace('\n', ''), '\n', old_bondtypes_line)
									cont=True
							if bondtype_line.split()[0] == old_bondtypes_line.split()[1] and bondtype_line.split()[1] == old_bondtypes_line.split()[0]:
								if "{:.3f}".format(float(bondtype_line.split()[3])) == "{:.3f}".format(float(old_bondtypes_line.split()[3])) and "{:.3f}".format(float(bondtype_line.split()[4])) == "{:.3f}".format(float(old_bondtypes_line.split()[4])):
									cont=True
								else:
									print('not copied:   \n',bondtype_line.replace('\n', ''), '\n', old_bondtypes_line)
									cont=True
				if cont==False:
					new_bonded.write(bondtype_line)
				cont=False
		if angle_types==True and len(bonded_line.replace('\n', ''))==0:
			for angletype_line in angletypes:
				if any(angletype_line.split()[0] in s for s in old_angle) or any(angletype_line.split()[1] in s for s in old_angle) or any(angletype_line.split()[2] in s for s in old_angle):
					for old_angle_line in old_angle:
						if old_angle_line[0][0] not in (';', '#', '[',''):
							if angletype_line.split()[0] == old_angle_line.split()[0] and angletype_line.split()[1] == old_angle_line.split()[1] and angletype_line.split()[2] == old_angle_line.split()[2]:
								if "{:.3f}".format(float(angletype_line.split()[4])) == "{:.3f}".format(float(old_angle_line.split()[4])) and "{:.3f}".format(float(angletype_line.split()[5])) == "{:.3f}".format(float(old_angle_line.split()[5])):
									cont=True
								else:
									print('not copied:   \n',angletype_line.replace('\n', ''), '\n', old_angle_line)
									cont=True
							if angletype_line.split()[2] == old_angle_line.split()[2] and angletype_line.split()[1] == old_angle_line.split()[1] and angletype_line.split()[0] == old_angle_line.split()[0]:
								if "{:.3f}".format(float(angletype_line.split()[4])) == "{:.3f}".format(float(old_angle_line.split()[4])) and "{:.3f}".format(float(angletype_line.split()[5])) == "{:.3f}".format(float(old_angle_line.split()[5])):
									cont=True
								else:
									print('not copied:   \n',angletype_line.replace('\n', ''), '\n', old_angle_line)
									cont=True
				if cont==False:
					new_bonded.write(angletype_line)
				cont=False

		if dih_types==True and len(bonded_line.replace('\n', ''))==0 and idih==1:
			for dihtype_line in dihtypes:
				if any(dihtype_line.split()[0] in s for s in old_dih) or any(dihtype_line.split()[1] in s for s in old_dih) or any(dihtype_line.split()[2] in s for s in old_dih):
					for old_dih_line in old_dih:
						if old_dih_line[0][0] not in (';', '#', '[',''):
							if dihtype_line.split()[0] == old_dih_line.split()[0] and dihtype_line.split()[1] == old_dih_line.split()[1] and dihtype_line.split()[2] == old_dih_line.split()[2] \
							and dihtype_line.split()[3] == old_dih_line.split()[3] and dihtype_line.split()[7] == old_dih_line.split()[7]:
								if "{:.3f}".format(float(dihtype_line.split()[4])) == "{:.3f}".format(float(old_dih_line.split()[4])) and "{:.3f}".format(float(dihtype_line.split()[5])) == "{:.3f}".format(float(old_dih_line.split()[5])):
									cont=True
								else:
									print('not copied:   \n',dihtype_line.replace('\n', ''), '\n', old_dih_line)
									cont=True
				if cont==False:
					new_bonded.write(dihtype_line)
				cont=False

		if idih_types==True and len(bonded_line.replace('\n', ''))==0 and idih==2:
			for idihtype_line in idihtypes:
				if any(idihtype_line.split()[0] in s for s in old_idih) or any(idihtype_line.split()[1] in s for s in old_idih) or any(idihtype_line.split()[2] in s for s in old_idih):
					for old_idih_line in old_idih:
						if old_idih_line[0][0] not in (';', '#', '[',''):
							if idihtype_line.split()[0] == old_idih_line.split()[0] and idihtype_line.split()[1] == old_idih_line.split()[1] and idihtype_line.split()[2] == old_idih_line.split()[2] \
							and idihtype_line.split()[3] == old_idih_line.split()[3]:
								if "{:.3f}".format(float(idihtype_line.split()[4])) == "{:.3f}".format(float(old_idih_line.split()[4])) and "{:.3f}".format(float(idihtype_line.split()[5])) == "{:.3f}".format(float(old_idih_line.split()[5])):
									cont=True
								else:
									print('not copied:   \n',idihtype_line.replace('\n', ''), '\n', old_idih_line)
									cont=True
				if cont==False:
					new_bonded.write(idihtype_line)
				cont=False

parser = argparse.ArgumentParser()

parser.add_argument('-f', help='file to add vsites to',metavar='system.pdb',type=str)
parser.add_argument('-ligand', help='ligand to add molecule to',metavar='GLO',type=str)
parser.add_argument('-func', help='what to do, add to forcefild or add vsites to pdb',metavar='[forcefield, vsites]',type=str, choices= ['forcefield', 'vsites'])
parser.add_argument('-nter', help='add flag to make n terminus charged', action='store_true')
parser.add_argument('-cter', help='add flag to make c terminus charged', action='store_true')



args = parser.parse_args()

if args.func=='forcefield':
	os.system('mkdir new_forcefield')
	os.chdir('new_forcefield')
	os.system('cp -r /sansom/s136/bioc1534/Documents/scripts/forcefield_scripts/charmm36-jul2017-updated.ff .')
	os.system('cp ../'+args.ligand+'.itp .')
	os.system('cp ../charmm36.itp .')



	print('vsites might need to be added to merged.vsd file')
	print('constraint may need to be added to ffbonded.itp')

	ligand, charmm=args.itp+'.itp','charmm36.itp'
	forcefield='charmm36-jul2017-updated.ff/'

	moleculetype=input(ligand, 'moleculetype', False)
	atoms=input(ligand, 'atoms', False)
	bonds=input(ligand, 'bonds', False)
	idihtypes_lig=input(ligand, 'dihedrals', True)

	bondtypes=input(charmm, 'bondtypes', False)
	angletypes=input(charmm, 'angletypes', False)
	dihtypes=input(charmm, 'dihedraltypes', False)
	idihtypes=input(charmm, 'dihedraltypes', True)

	merged_add(forcefield, 'merged.rtp')
	nonbonded_add(forcefield, 'ffnonbonded.itp')
	bonded_add(forcefield, 'ffbonded.itp')


elif args.func=='vsites':
	os.system('mkdir add_vsites')
	os.chdir('add_vsites')
	os.system('cp -r /sansom/s136/bioc1534/Documents/scripts/forcefield_scripts/charmm36-jul2017-updated.ff .')
	os.system('cp ../'+args.f+' .')
	pdb2gmx()
	make_system()












