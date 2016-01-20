import os
#rep=("namespace ParaMEDMEM","namespace MEDCoupling")
#rep=("ParaMEDMEM::","MEDCoupling::")
#rep=("ParaMEDMEMImpl::","MEDCouplingImpl::")

#rep=("_ParaMEDMEM__","_MEDCoupling__")
#rep=("ParaMEDMEM_","MEDCoupling_")
#rep=("ParaMEDMEMData","MEDCouplingData")
dirs=["MEDCoupling","MEDCoupling/Test","MEDLoader","MEDLoader/Swig","MEDLoader/Test","MEDPartitioner","MEDPartitioner/Test","MEDPartitioner_Swig","RENUMBER","RENUMBER_Swig","INTERP_KERNELTest","ParaMEDMEM","ParaMEDLoader","ParaMEDMEMTest","ParaMEDMEM_Swig"]
dirname=dirs[-1]
i=0
for fi in os.listdir(dirname):
    fi2=os.path.join(dirname,fi)
    if not os.path.isfile(fi2):
        continue
    f=file(fi2) ; lines=f.readlines() ; del f
    lines2=[line.replace(*rep) for line in lines]
    if lines2!=lines:
        i+=1
        f=file(fi2,"w") ; f.writelines(lines2) ; f.flush()
    pass

print i
