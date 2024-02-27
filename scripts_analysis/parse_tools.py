from .general_parameters import *

def calculate_GC_capacitances(sizes,alldata=None):
    if alldata is None:
        alldata=main_parser(sizes)

    for icat,cation in enumerate(alldata.keys()):#['interface_ion','Na_ads']:
      for isize,size in enumerate(alldata[cation].keys()):

        pot_v_q=alldata[cation][size]['GC']['Pot_v_q']
        coeff_e_V,d=curve_fit(lin_fun,pot_v_q[:,1],pot_v_q[:,0])

        #Changing units of charge to muC/cm2
        pot_v_q[:,0]*=float(coverages[size])/Aperatm*1.6022*1e3
        coeff,d=curve_fit(lin_fun,pot_v_q[:,1],pot_v_q[:,0])

        alldata[cation][size]['GC']['Cap_muF/cm2']=coeff[0]
        alldata[cation][size]['GC']['Cap_e/V']=coeff_e_V[0]


def main_parser(sizes):
    results={}
    add_slab_reference(results)
    for cationdir0  in os.listdir(home):
#    for cationdir0  in ['K_on_Pt']:
        cationdir=os.path.join(home,cationdir0)
        if not os.path.isdir(cationdir): continue
        if '_on_' not in cationdir0: continue
        cation=cationdir0.split('_')[0]
        results[cation]={}
        for size in sizes:
            if size not in os.listdir(cationdir): continue
            results[cation][size] = {}
            sizedir=cationdir+'/'+size
            for method in ['mean_pot','GC']:
                if method not in os.listdir(sizedir): continue
                basepath = sizedir+'/'+method+'/'
                results[cation][size][method] = {}
                read_trajfiles(results,basepath,size,method,cation)

            for method in ['CE']:
                if method not in os.listdir(sizedir):
                    continue
                basepath = sizedir+'/'+method+'/'
                repetitions = os.listdir(basepath)
                results[cation][size][method] = {}
                for rep in repetitions:
                    repetition = np.product([int(i) for i in rep.split('_')[1:]])
                    reppath = basepath+rep+'/'
                    read_trajfiles(results,reppath,size,method,cation,repetition)
    calculate_GC_capacitances(size,results)


#    out=open('parsed_data.pkl','wb')
#    pkl.dump(results,out)
#    out.close()

    return results

def read_trajfiles(results,basepath,size,method,cation,repetition=None):
    states = os.listdir(basepath)
    res=results[cation][size][method]
    add_slab_reference(results)
    if method == 'CE':
            res[repetition]={}
    for state in states:
        statepath=basepath+state
        if not os.path.isdir(statepath): continue

        #if method != 'CE':
        #    res[state.split('_')[-1]]={}
        #else:
        if method == 'CE':
            res[repetition][state.split('_')[-1]]={}
        trajfiles = os.listdir(statepath)

        E_v_pot,pots,dipmoms=[],[],[]
        for traj in trajfiles:
            if (traj.split('.')[-1] != 'traj' or
                traj.split('_')[0] != 'sp'):
                continue
            atoms=read(statepath+'/'+traj)
            ne = atoms.calc.results['ne']
            pot= atoms.calc.results['electrode_potential']
            E=atoms.get_potential_energy()
            dip=atoms.get_dipole_moment()[2]
            print(statepath,method)
            #if method == 'mean_pot' and abs(ne) > 1e-4: continue
            print(statepath,ne)

            pots.append([-ne,pot])
            dipmoms.append([-ne,atoms.get_dipole_moment()[2]])
            E_v_pot.append([pot,-ne,E])
#            if method == 'mean_pot':
#                Es.append([-ne,E])
#                res[state.split('_')[-1]][-ne]=np.array([E,pot,dip])
#            elif method == 'GC':
#                Es.append([pot,E])
                #Charge not electrons are written into the database
#                res[state.split('_')[-1]][pot]=np.array([E,-ne])
#            else:
#                Es.append([pot,E])
#                res[repetition][state.split('_')[-1]]=np.array([pot,E])
        if method in ['mean_pot','GC']:
            res['Pot_v_q'] = np.array(pots)
            res['dipmom_v_q'] = np.array(dipmoms)
            res['E_v_pot_and_q'] = np.array(E_v_pot)
        else:
            res[repetition][state.split('_')[-1]]['E_v_pot_and_q'] = np.array(E_v_pot)


def add_slab_reference(res,reference_path='Pt_slab_only'):
    res['slab']={}
    for size in sizes:
        sizedir=home+'/'+reference_path+'/'+size
        if not os.path.isdir(sizedir): continue
        siz=res['slab'][size]={}
        for method in os.listdir(sizedir):
            methoddir=sizedir+'/'+method
            if not os.path.isdir(methoddir): continue
            trajfiles = os.listdir(methoddir+'/state_slab/')
            met=siz[method]={}
            Es,pots,dipmoms=[],[],[]
            for traj in trajfiles:
                if (traj.split('.')[-1] != 'traj' or
                    traj.split('_')[0] != 'sp'):
                    continue
                trajpath = methoddir+'/state_slab/'+traj
                atoms=read(trajpath)
                ne = atoms.calc.results['ne']
                pot= atoms.calc.results['electrode_potential']
                E=atoms.get_potential_energy()
                dip=atoms.get_dipole_moment()[2]
                pots.append([-ne,pot])
                dipmoms.append([-ne,atoms.get_dipole_moment()[2]])
                Es.append([pot,-ne,E])
#                if method == 'mean_pot':
#                    Es.append([-ne,E])
#                elif method == 'GC':
#                    Es.append([pot,E])

                    #if method == 'mean_pot':
                    #    met[-ne]=[E,pot,dip]
                    #elif method == 'GC':
                    #    met[pot]=[E,-ne]

            if method in ['mean_pot','GC']:
                met['Pot_v_q'] = np.array(pots)
                met['dipmom_v_q'] = np.array(dipmoms)
                met['E_v_pot_and_q'] = np.array(Es)



def get_sizes():
    sizes=[]
    for size in os.listdir(home):
        if ( not os.path.isdir(home+'/'+size) or
            size == 'first_relax'):
            continue
        sizes.append(size)
    return sizes


    #Do meanpot analysis

def collect_results(sizes, alldata=None):
    if alldata is  None:
        alldata=main_parser(sizes)
    allcedpot_results={}

    for icat,cation in enumerate(alldata.keys()):
     allcedpot_results[cation]={}
     if cation == 'slab':continue

     for isize,size in enumerate(alldata[cation].keys()):
        mpresults=[]
        gcresults=[]

        #Calculate mean potential results
        mpdata=alldata[cation][size]['mean_pot']
        slabdat=alldata['slab'][size]['mean_pot']['E_v_pot_and_q']
        for icFS,chFS in enumerate(mpdata['E_v_pot_and_q'][:,1]):
#           if abs(chFS) > 1e-5: continue
           for icsl,chsl in enumerate(slabdat[:,1]):
               if np.around(chsl-chFS,3): continue
               mean_pot = (slabdat[icsl,0]+
                           mpdata['E_v_pot_and_q'][icFS,0]) / 2
               dE = (mpdata['E_v_pot_and_q'][icFS,2]-
                     slabdat[icsl,2]-refs[cation]+(mean_pot-redoxpots[cation]))
               if cation == 'Mg':
                   dE+=mean_pot-redoxpots[cation]
               mpresults.append([mean_pot,dE]) #,mpdata['FS'][chsl][2]])
        mpdata['dE_v_pot'] = mpresults

        #Calculate GC results

        gcdata=alldata[cation][size]['GC']
        isort=np.argsort(gcdata['E_v_pot_and_q'][:,0])
        gcdat=gcdata['E_v_pot_and_q'][isort]
        slabdat=alldata['slab'][size]['GC']['E_v_pot_and_q']
        dn,dE_canon=[],[]
        for ipFS,FSpot in enumerate(gcdat[:,0]):
            for ipsl,slpot in enumerate(slabdat[:,0]):
                if np.around(FSpot-slpot,1): continue
                if FSpot < 3.5: continue
                #dE_base=gcdat[ipFS,2]-slabdat[ipsl,2] - refs[cation]-redoxpots[cation]
                dE = gcdat[ipFS,2]-\
                      slabdat[ipsl,2]-refs[cation]+(FSpot-redoxpots[cation])

                #The subtraction on the right is opposite because of cahrge vs number of electrons
                dn.append([FSpot,1+(slabdat[ipsl,1]-gcdat[ipFS,1])])
                dE_canon.append([FSpot,dE - (1+(slabdat[ipsl,1]-gcdat[ipFS,1]))*FSpot])

                if cation == 'Mg':
                    dE+=FSpot-redoxpots[cation]
                gcresults.append([FSpot,dE])

        gcdata['dE_v_pot']=gcresults[1:]
        gcdata['dE_v_pot']=gcresults
        gcdata['dn_v_pot']=np.array(dn)
        gcdata['dE_canon_v_pot']=np.array(dE_canon)

        #assert (np.around(gcdata['dE_v_pot'][0][1] - gcdata['dE_canon_v_pot'][1,1] + gcdata['dn_v_pot'][1,1]*gcdata['dn_v_pot'][1,0],3))

        #Calculate extrapolation results
        if 'CE' not in alldata[cation][size].keys():
            print('CE data of %s and coverage %s is missing'%(cation,coverages[size]))
            continue
        cedata=alldata[cation][size]['CE']
        ceresults=[]
        for rep in cedata.keys():
            if ('FS' not in cedata[rep].keys() or
                'slab' not in cedata[rep].keys()): continue
            if not len(cedata[rep]['slab']['E_v_pot_and_q']) or \
                not len(cedata[rep]['FS']['E_v_pot_and_q']): continue

            pot=cedata[rep]['FS']['E_v_pot_and_q'][0,0]
            E_sl=cedata[rep]['slab']['E_v_pot_and_q'][0,2]
            E_FS=cedata[rep]['FS']['E_v_pot_and_q'][0,2]

            dE = E_FS-E_sl-refs[cation]+pot-redoxpots[cation]

            if cation == 'Mg':
                    dE+=pot-redoxpots[cation]

            dpot = pot-cedata[rep]['slab']['E_v_pot_and_q'][0,0]
            ceresults.append([dpot,dE])

        if not ceresults: continue
        ceresults=np.array(ceresults)
        if len(ceresults) > 1:
            coeff,dummy=curve_fit(lin_fun,ceresults[:,0],ceresults[:,1])
            alldata[cation][size]['CE']['dE_v_pot']=[pot,coeff[1]]
            allcedpot_results[cation][size]=ceresults

    slabdips={}
    for size in alldata['slab']:
        mpdat=alldata['slab'][size]['mean_pot']['dipmom_v_q']
        for ic,charge in enumerate(mpdat):
            if abs(charge[0]) < 1e-4:
                slabdips[size]=mpdat[ic][1]

    for icat, cation in enumerate(alldata):
        dat=[]
        for size in alldata[cation]:
            mpdat=alldata[cation][size]['mean_pot']['dipmom_v_q']
            for ic,charge in enumerate(mpdat):
                if abs(charge[0]) < 1e-4:
                    #dat.append([float(coverages[size]),mpdat[ic][1]-slabdips[size]])
                    alldata[cation][size]['mean_pot']['dmu']=mpdat[ic][1]-slabdips[size]
                    break
    out=open('parsed_data.pkl','wb')
    pkl.dump((alldata,allcedpot_results),out)
    out.close()
    return alldata,allcedpot_results


if __name__ == '__main__':
    main()
