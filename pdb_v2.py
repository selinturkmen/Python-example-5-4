def protein_parser(pdbfile):
 

    pdbfile =open("1ubq.pdb","r")
    atoms = []
    position = []
    
    for line in pdbfile:
        if line.startswith('ATOM'):
            line = line.split()
            element = line[11]
            dim = line[6:9]
            dimne = []
            residue = line[3]
            for x in dim:
                dimne.append(float(x))
            elements=[]
            residues = []
            residues.append(residue)
            elements.append(element)
            position = dimne + elements+ residues
            atoms.append(position)

    return atoms

def center_of_mass(atoms):
    
    mass={'C':12.01, 'O':16.00, 'H':1.008, 'N': 14.01, 'S':32.06}
    residues=("ALA","ARG","ASN","ASP","ASX","CYS","GLU","GLN","GLX","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL")


    mass_cm = 0
    mass_z = 0
    mass_x = 0
    mass_y = 0 
  
    for x in atoms:
           ri_x = x[0]
           ri_y = x[1]
           ri_z = x[2]
           mz = x[3]
           mi = mass[mz]
           mass_cm = mass_cm + mi
           mass_x = mass_x + ri_x * mi
           mass_y = mass_y + ri_y * mi
           mass_z = mass_z+ ri_z * mi
           cm_x = mass_x / mass_cm
           cm_y = mass_y / mass_cm
           cm_z = mass_z / mass_cm
           cm = [cm_x, cm_y, cm_z]
    return cm


def shift(atoms,vec):
   

    atomsnew = []
    for a in atoms:
        xdimension = a[0] + vec[0]
        ydimension = a[1] + vec[1]
        zdimension = a[2] + vec[2]
        isa = a[3]
        new = [xdimension,ydimension,zdimension,isa]
        atomsnew.append(new)
    return atomsnew


def writepdb(atomcoor,infile,outfile):
 
atoms=protein_parser('1ubq.pdb')


com=center_of_mass(atoms)
print(f"Center of mass of the protein is {com[0]:6.3f}, {com[1]:6.3f}, {com[2]:6.3f}")

for i in range(3):
    com[i] *= -1
atomsnew=shift(atoms,com)


comnew=center_of_mass(atomsnew[:])

print(f"New center of mass of the protein is {comnew[0]:6.3f}, {comnew[1]:6.3f}, {comnew[2]:6.3f}")


writepdb(atomsnew,infile="1ubq.pdb",outfile="1ubq_new.pdb")


for i in range(3):
    com[i] *= -1

atoms=protein_parser('1ubq_new.pdb')
atomsnew=shift(atoms,com)


com=center_of_mass(atoms)


print(f"After shifting back center of mass of the protein is {com[0]:6.3f}, {com[1]:6.3f}, {com[2]:6.3f}")


print("\n"*2)

com=center_of_mass(atoms)
print(f"Center of mass of ALA residues is {com[0]:6.3f}, {com[1]:6.3f}, {com[2]:6.3f}")

com=center_of_mass(atoms)
print(f"Center of mass of ALA & GLY residues is {com[0]:6.3f}, {com[1]:6.3f}, {com[2]:6.3f}")

com=center_of_mass(atoms)
print(f"Center of mass of C atoms is {com[0]:6.3f}, {com[1]:6.3f}, {com[2]:6.3f}")

com=center_of_mass(atoms)
print(f"Center of mass of all atoms is {com[0]:6.3f}, {com[1]:6.3f}, {com[2]:6.3f}")

com=center_of_mass(atoms)
print(f"Center of mass of C & N atoms is {com[0]:6.3f}, {com[1]:6.3f}, {com[2]:6.3f}")

com=center_of_mass(atoms)
print(f"Center of mass of C & N atoms of ALA & GLY residues is {com[0]:6.3f}, {com[1]:6.3f}, {com[2]:6.3f}")

