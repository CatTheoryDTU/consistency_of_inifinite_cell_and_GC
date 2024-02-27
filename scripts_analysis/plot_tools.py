from .general_parameters import *

def create_dE_v_pot_plot(alldata,name='figure4.pdf',cations=['H','Li','Na','K','Cs'],xlim=[-0.8,+1.2],plot_r2=True,
        outdir='results/dE_vs_pot/',SHE_potential=4.4,potscale='SHE'):

    figall,axall=plt.subplots()
    figsingle,axsingle=plt.subplots(nrows=2,ncols=3,figsize=(12,8))

    #For legend
    axall.plot(np.nan,np.nan,'o',label='Mean potential',markeredgecolor='k',markerfacecolor='w')
    axall.plot(np.nan,np.nan,'d',label='Grand canonical',markeredgecolor='k',markerfacecolor='w')
    axall.plot(np.nan,np.nan,'^',label='Infinite cell size',markeredgecolor='k',markerfacecolor='w')
    if len(cations)%2:
        axsingle[-1][-1].plot(np.nan,np.nan,'o',label='Mean potential',markeredgecolor='k',markerfacecolor='w')
        axsingle[-1][-1].plot(np.nan,np.nan,'d',label='Grand canonical',markeredgecolor='k',markerfacecolor='w')
        axsingle[-1][-1].plot(np.nan,np.nan,'^',label='Infinite cell size',markeredgecolor='k',markerfacecolor='w')
    else:
        axsingle[0][-1].plot(np.nan,np.nan,'o',label='Mean potential',markeredgecolor='k',markerfacecolor='w')
        axsingle[0][-1].plot(np.nan,np.nan,'d',label='Grand canonical',markeredgecolor='k',markerfacecolor='w')
        axsingle[0][-1].plot(np.nan,np.nan,'^',label='Infinite cell size',markeredgecolor='k',markerfacecolor='w')
    xpan=0

    for icat,cation in enumerate(cations):
     figalone,axalone=plt.subplots()
     axalone.plot(np.nan,np.nan,'o',label='Mean potential',markeredgecolor='k',markerfacecolor='w')
     axalone.plot(np.nan,np.nan,'d',label='Grand canonical',markeredgecolor='k',markerfacecolor='w')
     axalone.plot(np.nan,np.nan,'^',label='Infinite cell size',markeredgecolor='k',markerfacecolor='w')
     thiscatsizes=sizes
     if cation in ['K']:
         thiscatsizes=sizes[1:]
     if cation in ['Cs']:
         thiscatsizes=sizes[2:]
     fullfit_data=[]

     if not icat%3  and icat:
         xpan+=1

     axsing=axsingle[xpan][icat%3]

     for isize,size in enumerate(reversed(thiscatsizes)):
        axalone.plot(np.nan,np.nan,'s',color=colors[isize],label=r'$\theta$='+coverages_for_legends[size])
        if cation == 'Na':
            axall.plot(np.nan,np.nan,'s',color=colors[isize],label=r'$\theta$='+coverages_for_legends[size])
            if len(cations)%2:
                axsingle[-1,-1].plot(np.nan,np.nan,'s',color=colors[isize],label=r'$\theta$='+coverages_for_legends[size]+' ML')
            else:
                axsingle[0,-1].plot(np.nan,np.nan,'s',color=colors[isize],label=r'$\theta$='+coverages_for_legends[size]+' ML')

        siz=alldata[cation][size]
        mark=['o','d','^']
        for im,method in enumerate(['mean_pot','GC','CE']):
            if method not in siz.keys(): continue
            if 'dE_v_pot' not in siz[method]: continue
            if not len(siz[method]['dE_v_pot']): continue
            if 1:
                plotdat=np.array(siz[method]['dE_v_pot'])

                if method == 'CE':
                    x,y=plotdat[0],plotdat[1]
                    fullfit_data.append(siz[method]['dE_v_pot'])
                else:
                    for i in siz[method]['dE_v_pot']:
                        fullfit_data.append([i[0],i[1]])
                    x,y=plotdat[:,0],plotdat[:,1]
                if potscale == 'SHE': x-=SHE_potential
                for thisax in [axall,axsing,axalone]:
                        thisax.plot(x,y,mark[im],
                            color=colors[isize],linewidth=0.5,markeredgecolor='k',markersize=markersize)
            else:
                print(method+' of %s and size %s is missing'%(cation,size))

        axsing.set_xticks([3.5,4,4.5,5,5.5])
        if potscale == 'SHE':
            axsing.set_xticks([-1.0,-0.5,0,0.5,1])


        if cation in ['H']:
            axsing.set_yticks(np.arange(-2,5,0.6))
        elif cation in ['Li']:
            axsing.set_yticks(np.arange(-2,5,0.5))
        elif cation in ['Na']:
            axsing.set_yticks(np.arange(-2,5,0.3))
        elif cation in ['K']:
            axsing.set_yticks(np.arange(-2,5,0.25))
        else:
            axsing.set_yticks(np.arange(-2,5,0.2))

     fullfit_data=np.array(fullfit_data)
     if potscale == 'SHE': fullfit_data[:,0]-=SHE_potential
     coeff,pcov=curve_fit(lin_fun,fullfit_data[:,0],fullfit_data[:,1])
     print(cation+':')
     print('Fitting curve coefficients (A*x+B):', coeff)
     print('Plating potential (absolute,vs SHE):', -coeff[1]/coeff[0], -coeff[1]/coeff[0]-SHE_potential)
     print('-'*13)

     texty=coeff[0]*0.5+coeff[1]-0.2
     cationy=coeff[0]*5.5+coeff[1]-0.1

     for thisax in [axall,axsing,axalone]:
         thisax.plot(xlim,coeff[0]*np.array(xlim)+coeff[1],'--',color='k')

     for thisax in [axall,axsing,axalone]:
         thisax.annotate(r'$\frac{d\Delta G}{d\Phi}$=%1.2fe'%coeff[0],(-0.5,texty),rotation=43,fontsize=24).draggable()

     plt.xlim(xlim)

     textsqy=coeff[0]*3.5+coeff[1]
     if plot_r2:
        residuals = fullfit_data[:,1] - lin_fun(fullfit_data[:,0], *coeff)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((fullfit_data[:,1]-np.mean(fullfit_data[:,1]))**2)
        r_squared = 1 - (ss_res / ss_tot)
        for thisax in [axall,axalone,axsing]:
            thisax.annotate(r'R$^2$=%1.2f'%r_squared,(0.3,min(fullfit_data[:,1])-0.1),fontsize=20,va='top').draggable()#=True)


        #axsing.annotate(r'R$^2$=%1.2f'%r_squared,(4.7,textsqy),fontsize=16).draggable()#=True)

#     axsing.annotate(r'$\Delta$ E=%1.2fe$\Phi$+%1.2f'%(coeff[0],coeff[1]),(4.0,textsqy+0.2),fontsize=16).draggable()
#     axalone.annotate(r'$\Delta$ E=%1.2fe$\Phi$+%1.2f'%(coeff[0],coeff[1]),(4.0,textsqy+0.2),fontsize=16).draggable()
     axalone.set_xlabel('(Mean) Potential $\phi$ [V vs SHE]')
     axalone.set_ylabel('$\Delta$G$_\phi$ [eV]')
     figalone.legend(fontsize=16,bbox_to_anchor=(0.95,0.405),ncol=2)
     figalone.tight_layout()
     figalone.savefig(outdir+'%s.pdf'%cation)
     #plt.show()
     figalone.gca()

    if len(cations)%2:
        axsingle[-1][-1].set_axis_off()
        axsingle[-1][-1].legend(loc='upper center',bbox_to_anchor=(0.5,0.9),fontsize=16)#,bbox_to_anchor=(0,0))
        #axsingle[-1][-1].legend(loc='upper center',bbox_to_anchor=(-2,1.05,3,0.3),fontsize=16)#,bbox_to_anchor=(0,0))
    else:
        axsingle[0][-1].legend(fontsize=14,bbox_to_anchor=(1,1))
    figall.legend(fontsize=14,bbox_to_anchor=(0.95,0.4),ncol=2)
    axall.set_xlabel('(Mean) Potential [V vs SHE]')
    axall.set_ylabel('Adsorption energy $\Delta G$ [eV]')

#    figsingle.legend(loc='upper center',bbox_to_anchor=(1,1,1,0.3),fontsize=16)
    axsingle[-1, 0].set_xlabel('.', color=(0, 0, 0, 0))
    axsingle[-1, 0].set_ylabel('.', color=(0, 0, 0, 0))
#    refpos=axsingle[1][0].get_position()
#    axsingle[1][0].set_position([0.5,0.15,refpos.width,refpos.height])
#    axsingle[1][1].set_position([0.6,0.15,refpos.width,refpos.height])
    figsingle.text(0.25,0.03,'(Mean) Potential $\phi$ / V vs SHE',ha='left',fontsize=20)
    figsingle.text(0.7,0.5,'(Mean) Potential $\phi$ / V vs SHE',ha='left',fontsize=20)
    figsingle.text(0.01,0.5,'$\Delta$G / eV',rotation='vertical',va='center',fontsize=25)
    cationy=coeff[0]*5.5+coeff[1]-0.1
    figsingle.text(0.1,0.95,'*H',va='center',fontsize=30,fontweight='bold')
    figsingle.text(0.415,0.95,'*Li',va='center',fontsize=30,fontweight='bold')
    figsingle.text(0.73,0.95,'*Na',va='center',fontsize=30,fontweight='bold')
    figsingle.text(0.1,0.485,'*K',va='center',fontsize=30,fontweight='bold')
    figsingle.text(0.415,0.485,'*Cs',va='center',fontsize=30,fontweight='bold')
    figall.tight_layout()
    figsingle.tight_layout()
    figsingle.savefig(outdir+'separate_panels.pdf')
    figall.savefig(outdir+'singlepanel.pdf')
#    figsingle.show()
    #plt.savefig(name)
    plt.show()
    #plt.close()



def plot_ce_dpot_plot(ceresults,allresults,cations=['Li','Na','K','Cs']):
  fontsize=20
  for icat,cation in enumerate(cations):
    fig,ax = plt.subplots(1,2,sharey=True,figsize=(15,5),gridspec_kw={'width_ratios': [3, 1]})

    de_v_pot=[]
    for isize,size in enumerate(reversed(ceresults[cation])):
        res=ceresults[cation][size]
        sorti=np.argsort(res[:,0])
        coeff,dummy=curve_fit(lin_fun,res[:,0],res[:,1])
        xvals= [i for i in res[:,0]]+[0.0]
        ax[0].plot(xvals,coeff[0]*np.array(xvals)+coeff[1],'-',color=colors[isize])
        ax[0].plot(res[sorti][:,0],res[sorti][:,1],'o',color=colors[isize],markeredgecolor='k',markersize=markersize)
        de_v_pot.append([allresults[cation][size]['CE']['dE_v_pot'][0],coeff[1]])

        if cation not in ['K','Cs']:
         if size in allresults[cation].keys():
            if size=='sqrt13xsqrt13':
             ax[0].annotate(r'$\frac{d\Delta G}{d\Delta \Phi}$ = %1.2fe'%coeff[0], (xvals[0]-0.08,coeff[0]*(xvals[0]+0.0)+coeff[1]),rotation=coeff[0]*70,fontsize=fontsize*0.75).draggable()
             ax[0].annotate(r'$\theta$=%s,$\phi$=%1.2f'
                     %(coverages_for_legends[size],allresults[cation][size]['CE']['dE_v_pot'][0]),
                     (0.01,coeff[1]-0.04),fontsize=fontsize).draggable()
    #            ax[0].annotate(r'$\frac{d\Delta G}{d\Delta \Phi}$ = %1.2fe'%coeff[0], (xvals[0]-0.02,coeff[0]*(xvals[0]-0.14)+coeff[1]),rotation=coeff[0]*70).draggable()
            else:
             ax[0].annotate(r'$\frac{d\Delta G}{d\Delta \Phi}$ = %1.2fe'%coeff[0], (xvals[0]-0.05,coeff[0]*(xvals[0])+coeff[1]),rotation=coeff[0]*70,fontsize=fontsize*0.75).draggable()
             ax[0].annotate(r'$\theta$=%s,$\phi$=%1.2f'%(coverages_for_legends[size],
                 allresults[cation][size]['CE']['dE_v_pot'][0]),
                 (0.01,coeff[1]-0.02),fontsize=fontsize).draggable()

    de_v_pot=np.array(de_v_pot)
    coeff,dummy=curve_fit(lin_fun,de_v_pot[:,0],de_v_pot[:,1])
    ax[1].plot(de_v_pot[:,0],coeff[0]*de_v_pot[:,0]+coeff[1],'-k')
    for isize,size in enumerate(reversed(ceresults[cation])):
        ax[1].plot(de_v_pot[isize,0],de_v_pot[isize,1],'o',color=colors[isize],markeredgecolor='k',markersize=markersize)
    x_ano=de_v_pot[-1,0]+0.1
    y_ano=coeff[0]*x_ano+coeff[1]+0.05
    ax[1].annotate(r'$\frac{d\Delta G^\infty_n}{d\Phi}$ = %1.2fe'%coeff[0],(x_ano,y_ano),rotation=coeff[0]*65,fontsize=20)

    ax[0].axvline(x=0,linestyle=':',color='k')
    ax[0].set_xlim([-1.1,0.24])
    ax[0].set_xlabel('Potential mismatch $\Delta \phi$ [V]')
    ax[1].set_xlabel('Potential $\phi$ [V]')
    ax[0].set_ylabel('$\Delta$G$_n$ [eV]')
    fig.tight_layout()
#    plt.savefig('results/'+cation+'_CE_E_vs_dpot.pdf')
    plt.show()
    plt.close()

def plot_GC_dEcanondphi_vs_cap(alldata,cations=['slab','H','Li','Na','K','Cs'],outdir='results/',show=False):

    outnames=['dG_can_dphi_in_GC_v_Cap.pdf']
    varnames=['E_v_pot_and_q']
    ylab=[r'$\frac{\partial G_{\mathrm{A*}/*}}{\partial \phi}[\phi]$ [eV]']

    allslopes=[]
    for ivar,var in enumerate(varnames):
        for icat, cation in enumerate(cations):
            plt.plot(np.nan,np.nan,color=catcolors[icat], marker='s',label=cation)
            for isize,size in enumerate(alldata[cation]):
                if icat==len(cations)-1:
                    plt.plot(np.nan,np.nan,color='w',marker=markers[isize],label=coverages_for_legends[size]+'ML',markeredgecolor='k')
                y=np.array(alldata[cation][size]['GC'][var][:,2])+np.array(alldata[cation][size]['GC'][var][:,1])*np.array(alldata[cation][size]['GC'][var][:,0])
                x=np.array(alldata[cation][size]['GC'][var][:,0])
                cap=-alldata[cation][size]['GC']['Cap_e/V']
                linear=1
                if linear:
                    coeff,d=curve_fit(lin_fun,x,y)
                # minus sign before y-vals come from the fact that we want electrons not charges
                    plt.plot(cap,coeff[0],marker=markers[icat],color=colors[isize],markeredgecolor='k')
                    allslopes.append([cap,coeff[0]])
                else:

                    coeff,d=curve_fit(quad_fun,x,y)
                # minus sign before y-vals come from the fact that we want electrons not charges
                    seconderiv=False
                    for pot in np.arange(4,5,5):
                        if seonderiv:
                            y=2*coeff[0]
                        else:
                            y=2*coeff[0]*pot+coeff[1]
                        plt.plot(cap,2*coeff[0],marker=markers[icat],color=colors[isize],markeredgecolor='k')
                    allslopes.append([cap,2*coeff[0]*pot+coeff[1]])

        allslopes=np.array(allslopes)
        coeff,d=curve_fit(lin_fun,allslopes[:,0],allslopes[:,1])
        plt.plot([min(allslopes[:,0]),max(allslopes[:,0])],[min(allslopes[:,0])*coeff[0]+coeff[1],max(allslopes[:,0])*coeff[0]+coeff[1]],'--k',label='Slope: %1.2f'%coeff[0])
        plt.ylabel(ylab[ivar])
        plt.xlabel('C$_{\mathrm{A}*/*}$ [e$^-$/V]')
        plt.legend()
        plt.tight_layout()
        if show: plt.show()
        else: plt.savefig(outdir+outnames[ivar])
        plt.close()

def plot_GC_dne_vs_pot(alldata,cations=['H','Li','Na','K','Cs'],outdir='results/',show=False):

    outnames=['GC_dn_vs_pot.pdf','dG_can_in_GC_v_pot.pdf','dG_GC.pdf','G_can_in_GC_v_pot.pdf']
    varnames=['dn_v_pot','dE_canon_v_pot','dE_v_pot','E_v_pot_and_q']
    ylab=['(1+$\Delta$n) [e$^-$]',
            'G$_{\mathrm{A*}}[\phi]$-G$_{\mathrm{*}}[\phi]-\mu^0_\mathrm{A}$ [eV]',
            '$\Delta$G[$\phi$] [eV]','G$_{\mathrm{A*}}[\phi]$']

    for ivar,var in enumerate(varnames):
        for icat, cation in enumerate(cations):
            plt.plot(np.nan,np.nan,'s',color=catcolors[icat],label=cation)
            for isize,size in enumerate(alldata[cation]):
                if icat==len(cations)-1:
                    plt.plot(np.nan,np.nan,marker=markers[isize],color='w',markeredgecolor='k',label=coverages_for_legends[size]+'ML')
                dat=np.array(alldata[cation][size]['GC'][var])
                dat=dat[np.argsort(dat[:,0])]
                # minus sign before y-vals come from the fact that we want electrons not charges
                plt.plot(dat[:,0],dat[:,1],'-.',marker=markers[isize],color=catcolors[icat],markeredgecolor='k')
        plt.ylabel(ylab[ivar])
        plt.xlabel('Potential $\phi$ [V]')
        plt.legend()
        plt.tight_layout()
        if show: plt.show()
        else: plt.savefig(outdir+outnames[ivar])
        plt.close()


def plot_dipoles_from_mean_potential(alldata,cations=['H','Li','Na','K','Cs'],outfile='results/dipole_moments_vs_coverage.pdf',show=False):

    pltdat={}

    for icat, cation in enumerate(cations):
        pltdat[cation]=[]
        for size in alldata[cation]:
            pltdat[cation].append([coverages[size],alldata[cation][size]['mean_pot']['dmu']])

        pltdat[cation] = np.array(pltdat[cation])
        plt.plot(pltdat[cation][:,0],pltdat[cation][:,1],'o--',label=cation,markeredgecolor='k')
    plt.ylabel('Surface dipole per adsorbate [e$\mathrm{\AA{}}$]',fontsize=20)
    plt.xlabel('Cation coverage')
    plt.legend()
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

    for icat, cation in enumerate(cations):
        pltdat[cation]=[]
        for size in alldata[cation]:
            dn=np.array(alldata[cation][size]['GC']['dn_v_pot'])
            for i in range(len(dn)):
                pltdat[cation].append([alldata[cation][size]['mean_pot']['dmu'],dn[i,1]])

        pltdat[cation] = np.array(pltdat[cation])
        plt.plot(pltdat[cation][:,0],pltdat[cation][:,1],'o--',label=cation,markeredgecolor='k')
    plt.xlabel('Surface dipole per adsorbate [e$\mathrm{\AA{}}$]',fontsize=20)
    plt.ylabel('1+$\Delta n$ [e]')
    plt.legend()
    plt.tight_layout()
    if show: plt.show()
    else:
        plt.savefig('results/dn_vs_dipole_moments.pdf')
    plt.close()

def plot_slope_vs_elneg(alldata,cations=['H','Li','Na','K','Cs']):
#def plot_slope_vs_elneg(alldata,cations=['Li','Na','K','Cs']):
    pltdat={}
    figeleneg,axeleneg=plt.subplots()
    slopes=[]
    for icat,cation in enumerate(cations):
     thiscatsizes=list(alldata[cation].keys())
     fullfit_data=[]
     if cation in ['K']:
         thiscatsizes=thiscatsizes[1:]
     if cation in ['Cs']:
         thiscatsizes=thiscatsizes[2:]
     for isize,size in enumerate(reversed(thiscatsizes)):
        siz=alldata[cation][size]
        if len(siz['mean_pot']):
            fullfit_data.extend([i for i in siz['mean_pot']['dE_v_pot']])
        if len(siz['GC']):
            fullfit_data.extend([i for i in siz['GC']['dE_v_pot'] if i[0] > 3.5])
        if 'CE' not in siz.keys(): continue
#        print(cation,size,siz['CE'].keys(),fullfit_data)
#        dd
        if 'dE_v_pot' not in siz['CE']:
            print(f'Couldnt find CE results for {cation} {size}')
            continue
        fullfit_data.append([siz['CE']['dE_v_pot'][0],siz['CE']['dE_v_pot'][1]])# for i in siz['CE']])

     fullfit_data=np.array(fullfit_data)
     coeff,pcov=curve_fit(lin_fun,fullfit_data[:,0],fullfit_data[:,1])
     if cation != 'H':
         slopes.append([elneg[cation],coeff[0]])

     axeleneg.plot(elneg[cation],coeff[0],'sk',markeredgecolor='k',markersize=10)
     axeleneg.axhline(y=1,linestyle='--',color='k')
     if cation == 'Li':
         axeleneg.annotate(cation,(elneg[cation]-0.017,coeff[0]-0.005), fontsize=30)
     else:
         axeleneg.annotate(cation,(elneg[cation]+0.005,coeff[0]+0.005), fontsize=30)

    slopes=np.array(slopes)
    coeff,pcov=curve_fit(lin_fun,slopes[:,0],slopes[:,1])
#    axeleneg.plot(slopes[:,0],slopes[:,0]*coeff[0]+coeff[1],'--k')
#    axeleneg.annotate(r'y=%1.2fx%1.2f'%(coeff[0],coeff[1]),(0.93,slopes[-1][1]),fontsize=24)

    axeleneg.set_xlabel('Pauling electronegativity')
    axeleneg.set_ylabel(r'$\frac{d\Delta G}{d\Phi} [\mathrm{e}]$')
    figeleneg.tight_layout()
    figeleneg.savefig('results/Slope_vs_elneg.pdf')



def plot_GC_capacitances(outdir='results/Capacitances/',alldata=None,Aperatm=6.928203230275509):

    plt.close()

    if alldata is None:
        alldata=collect_results(sizes)

    figgcsingleppanel,axgcsinglepanel=plt.subplots()
    allcaps={}
    for icat,cation in enumerate(alldata.keys()):
      if cation == 'slab': continue
      figgc,axgc=plt.subplots()
      figmp,axmp=plt.subplots()
      fig_Cvsize,ax_Cvsize=plt.subplots()
      fig_Cvcov,ax_Cvcov=plt.subplots()
      mpdat=[]
      caps={}
      for isize,size in enumerate(alldata[cation].keys()):
        gcdata=alldata[cation][size]['GC']
        slabdata=alldata['slab'][size]['GC']
#        mpdata=alldata[cation][size]['mean_pot']
        markers={'slab':'s','FS':'o'}
        for state in ['slab','FS']:
            if state not in caps:
                caps[state]=[]

            pltdat=[]
            if state == 'FS':
                caps[state].append([float(coverages[size]),float(gcdata['Cap_muF/cm2'])])
                pltdat=gcdata['E_v_pot_and_q'][np.argsort(gcdata['E_v_pot_and_q'][:,0])][:,:2]
            else:
                caps[state].append([float(coverages[size]),float(slabdata['Cap_muF/cm2'])])
                pltdat=slabdata['E_v_pot_and_q'][np.argsort(gcdata['E_v_pot_and_q'][:,0])][:,:2]

            pltdat=np.array(pltdat)
            pltdat[:,1]*=float(coverages[size])/Aperatm*1.6022*1e3
            axgc.plot(pltdat[:,0],pltdat[:,1],markers[state],color=colors[isize])

        axgc.plot(np.nan,np.nan,'s',color=colors[isize],label=r'$\theta$=%s'%(coverages_for_legends[size]))

      #mpdat=np.array(mpdat)
      #axmp.plot(mpdat[:,0],mpdat[:,1],'o-')

      axgc.set_xlabel('Potential [V]')
      axgc.set_ylabel('Applied charge density [$\mu C/cm^2$]')
      axmp.set_xlabel('Potential [V]')
      axmp.set_ylabel('E [eV]')
      figgc.legend()
      figgc.tight_layout()
      figgc.savefig(outdir+'Charge_vs_potential_%s.pdf'%cation)
      figmp.savefig(outdir+'Energy_vs_potential_mean_pot_%s.pdf'%cation)
      plt.close()

      for state in caps:
          caps[state]=np.array(caps[state])
      ax_Cvsize.plot(1/caps['slab'][:,0],caps['slab'][:,1],'ks-',label='Bare slab')
      ax_Cvsize.plot(1./caps['FS'][:,0],caps['FS'][:,1],'ro-',label='*%s'%cation)
      ax_Cvsize.set_xlabel('Surface area [atoms]')
      ax_Cvsize.set_ylabel('Capacitance [muF/cm2]')
      fig_Cvsize.legend()
      fig_Cvsize.tight_layout()
      fig_Cvsize.savefig(outdir+'Capacitance_vs_size_%s.pdf'%cation)

      ax_Cvcov.plot(caps['slab'][:,0],caps['slab'][:,1],'ks-',label='Bare slab')
      ax_Cvcov.plot(caps['FS'][:,0],caps['FS'][:,1],'ro-',label='*%s'%cation)
      ax_Cvcov.set_xlabel('Coverage')
      ax_Cvcov.set_ylabel('Capacitance [muF/cm2]')
      fig_Cvcov.legend()
      fig_Cvcov.tight_layout()
      fig_Cvcov.savefig(outdir+'Capacitance_vs_coverage_%s.pdf'%cation)

      #plot capacitance parity plot
      plt.plot([15,30],[15,30],'--k')
      plt.plot(caps['slab'][:,1],caps['FS'][:,1],'ks-')
      plt.xlabel('Capacitance of * + (%s$^+_{aq}$ + e$^-$) [muF/cm$^2$]'%cation)
      plt.ylabel('Capacitance of *%s [muF/cm$^2$]'%cation)
      plt.tight_layout()
      plt.savefig(outdir+'Capacitance_parity_plot_%s.pdf'%cation)
      plt.close()

    plt.close()
    fig,ax=plt.subplots()
    for icat,cation in enumerate(['slab','H','Li','Na','K','Cs']):
      pltdat=[]
      for size in alldata[cation]:
        pltdat.append([1/float(coverages[size]),alldata[cation][size]['GC']['Cap_muF/cm2']])
      pltdat=np.array(pltdat)
      ax.plot(pltdat[:,0],pltdat[:,1],
                'o--',label=cation,color=catcolors[icat],markeredgecolor='k')
    ax.set_xlabel('Coverage')
    ax.set_ylabel('Capacitance [muF/cm2]')
    fig.legend()
    fig.tight_layout()
    fig.savefig(outdir+'Capacitance_vs_coverage_all_cations.pdf')
    plt.close()

    fig,ax=plt.subplots()
    for icat,cation in enumerate(['slab','H','Li','Na','K','Cs']):
      pltdat=[]
      for size in alldata[cation]:
        pltdat.append([1/float(coverages[size]),alldata[cation][size]['GC']['Cap_muF/cm2']])
      pltdat=np.array(pltdat)
      ax.plot(pltdat[:,0],pltdat[:,1],
                'o--',label=cation,color=catcolors[icat],markeredgecolor='k')
    ax.set_xlabel('Surface  area [atoms]')
    ax.set_ylabel('Capacitance [muF/cm2]')
    fig.legend()
    fig.tight_layout()
    fig.savefig(outdir+'Capacitance_vs_size_all_cations.pdf')
    plt.close()

def plot_dn_and_E_vs_Capmismatch(alldata,outdir='results/'):
    figdn,axdn=plt.subplots()
    figdndphi_phi,axdndphi_phi=plt.subplots()
    figdce,axdce=plt.subplots()
    figdg,axdg=plt.subplots()
    alldnslopes,alldceslopes,alldndphi_phislopes=[],[],[]
    dg_vs_analytical_slope=[]
    for icat,cation in enumerate(['H','Li','Na','K','Cs']):
        for isize,size in enumerate(alldata[cation]):
            for axs in [axdn,axdce,axdg,axdndphi_phi]:
                if not isize:
                    axs.plot(np.nan,np.nan,'s',color=catcolors[icat],label=cation)

                if cation == 'Cs':
                    axs.plot(np.nan,np.nan,markers[isize],markeredgecolor='k',
                            markerfacecolor='w',label=coverages_for_legends[size]+'ML')

            dn=alldata[cation][size]['GC']['dn_v_pot']
            dCE=alldata[cation][size]['GC']['dE_canon_v_pot']#+2
            dG=np.array(alldata[cation][size]['GC']['dE_v_pot'])

#            dCE[:,1]+=2*dn[:,1]
#            print(dCE[:,1]+dn[:,1]*dn[:,0],dG[:,1])
#            dd

            #Numerical dn/dphi for plotting dn/dphi*phi
            dndphi = np.diff(dn[:,1])/np.diff(dn[:,0])
            dndphi_phi = np.diff(dn[:,1])/np.diff(dn[:,0])*(dn[:-1,0]+np.diff(dn[:,0]))

            #Numerical dG_can/dphi for plotting dG_can2/dphi2, currently not plotted
            dpots=dCE[:-1,0]+np.diff(dCE[:,0])
            dGdphi = np.diff(dCE[:,1])/np.diff(dCE[:,0])
            dG2dphi2 = np.diff(dGdphi)/np.diff(dpots)*(dpots[:-1]+np.diff(dpots[:]))

            #dn/dphi from fit (assuming linear behaviour)
            coeff,dummy=curve_fit(lin_fun,dn[:,0],dn[:,1])
            print(dndphi_phi,coeff)
            #dd
            coeffdce,dummy=curve_fit(lin_fun,dCE[:,0],dCE[:,1])
            #coeffdce,dummy=curve_fit(quad_fun,dCE[:,0],dCE[:,1])
            coeffdg,dummy=curve_fit(quad_fun,dG[:,0],dG[:,1])

            plotx=-alldata[cation][size]['GC']['Cap_e/V']+alldata['slab'][size]['GC']['Cap_e/V']

            alldnslopes.append([plotx,coeff[0]])
            for i in range(len(dndphi_phi)):
                alldndphi_phislopes.append([plotx,dndphi_phi[i]])
            alldceslopes.append([plotx,coeffdce[0]])

            # Plot dn/dphi vs Capmismath, fig S4
            axdn.plot(plotx,coeff[0],
                        markers[isize]+'--',color=catcolors[icat],markeredgecolor='k')

            # Plot dG - dn/dphi*phi vs Capmismath, fig S4
            #axdndphi_phi.plot([plotx]*len(dndphi_phi),dndphi_phi+dGdphi,
            axdndphi_phi.plot(dndphi_phi,dGdphi,
                        markers[isize],color=catcolors[icat],markeredgecolor='k')
#            axdndphi_phi.plot([plotx]*(len(dn)-1),coeff[0]*(dn[:-1,0]+np.diff(dn[:,0])),#dndphi_phi,
#                        markers[isize]+'--',color=catcolors[icat],markeredgecolor='k')

            # Plot dG/dphi, fig S3
            for potline in dCE:
                pot=potline[0]
                if pot < 4.0 and cation == 'Cs': continue

                #axdce.plot(-pot*plotx,2*coeffdce[0]*pot+coeffdce[1]+coeff[0]*pot,
                axdce.plot(plotx,coeffdce[0],
                        markers[isize]+'--',color=catcolors[icat],markeredgecolor='k')

            #dg_vs_analytical_slope.append([1+coeffdce[0]+coeff[0]*4.4+dn[1,1],coeffdg[0]])
            # Plot dDG/dphi vs 1+ne, fig S5
            for potline in dn:
                pot=potline[0]
                axdg.plot(potline[1],(2*coeffdg[0]*pot+coeffdg[1]),
                    markers[isize],color=catcolors[icat],markeredgecolor='k')


    #Get slope of slopes:
    alldceslopes=np.array(alldceslopes)
    coeffdce,d=curve_fit(lin_fun,alldceslopes[:,0],alldceslopes[:,1])
    alldndphi_phislopes=np.array(alldndphi_phislopes)
    coeffdndphi_phi,d=curve_fit(lin_fun,alldndphi_phislopes[:,0],alldndphi_phislopes[:,1])

    axdn.plot([0,0.7],[0,0.7],':k')
    axdce.plot([0,0.7],[0,coeffdce[0]*0.7+coeffdce[1]],':k',label='Slope: %1.2f' %coeffdce[0])
#    axdndphi_phi.plot([0,0.7],[0,coeffdndphi_phi[0]*0.7+coeffdndphi_phi[1]],':k',label='Slope: %1.2f' %coeffdndphi_phi[0])
    axdndphi_phi.plot([0,3],[0,-3],':k')
    axdg.plot([-0.5,1],[-0.5,1],'--k')

    axdn.set_xlabel('C$_{*A}-C_{*}$ [e$^-$/V]')
    #axdndphi_phi.set_xlabel('C$_{*A}-C_{*}$ [e$^-$/V]')
    axdndphi_phi.set_xlabel(r'$\frac{\partial \Delta n}{\partial \phi}\phi$ [e$^-$]')
    axdndphi_phi.set_ylabel(r'$\frac{\partial (G_{*A}[\phi]-G_{*}[\phi]) }{\partial \phi}$ [e$^-$]')
    axdce.set_xlabel('C$_{*A}-C_{*}$ [e$^-$/V]')
    axdg.set_xlabel(r'1+$\Delta$n [e$^-$]')
    axdn.set_ylabel(r'$\frac{\partial \Delta n}{\partial \phi}$ [e$^-$/V]')
    axdce.set_ylabel(r'$\frac{\partial (G_{*A}[\phi]-G_{*}[\phi]) }{\partial \phi}$ [e$^-$]')
    axdg.set_ylabel(r'$\frac{\partial \Delta G[\phi]}{\partial \phi}$ [e$^-$]')

    figdn.legend(loc='lower left',bbox_to_anchor=(0.85,0.16))
    #figdndphi_phi.legend(loc='lower left',bbox_to_anchor=(0.85,0.16))
    figdndphi_phi.legend()
    figdce.legend()
    figdg.legend(loc='lower left',bbox_to_anchor=(0.8,0.2))
    figdn.tight_layout()
    figdndphi_phi.tight_layout()
    figdce.tight_layout()
    figdg.tight_layout()
    figdn.savefig(outdir+'dndphi_vs_Capmismatch_all_cations.pdf')
    figdndphi_phi.savefig(outdir+'dndphi_phi_vs_Capmismatch_all_cations.pdf')
    figdce.savefig(outdir+'dCanonical_energies_of_GC_vs_Capmismatch_all_cations.pdf')
    figdg.savefig(outdir+'dG_GC_versus_analytical_slope.pdf')
    plt.show()
    plt.close()

def plot_GC_dG_vs_Capratio(alldata,outdir='results/',seconderiv=True):
    fig,ax=plt.subplots()
    cations=['H','Li','Na','K','Cs']
    allslopes=[]
    for icat,cation in enumerate(cations):
        for isize,size in enumerate(alldata[cation]):
            if not isize:
                ax.plot(np.nan,np.nan,'s',color=catcolors[icat],label=cation)
            if cation == cations[-1]:
                ax.plot(np.nan,np.nan,markers[isize],markeredgecolor='k',
                        markerfacecolor='w',label=coverages_for_legends[size]+'ML')


            dn=np.array(alldata[cation][size]['GC']['dE_v_pot'])

            coeff,dummy=curve_fit(quad_fun,dn[:,0],dn[:,1])
            print(cation,size,coeff)
            plotx=- alldata[cation][size]['GC']['Cap_e/V'] + alldata['slab'][size]['GC']['Cap_e/V']
            seconderiv=True
            if seconderiv:
                slope=2*coeff[0]
                allslopes.append([plotx,slope])
                ax.plot(plotx,slope,
                        markers[isize]+'--',color=catcolors[icat],markeredgecolor='k')
            else:
                for potline in dn:
                    pot = potline[0]
                    slope=(2*coeff[0]*pot+coeff[1])#/alldata['slab'][size]['GC']['Cap_e/V']
                    allslopes.append([pot,slope])

                    ax.plot(pot,slope,
                        markers[isize]+'--',color=catcolors[icat],markeredgecolor='k')


    #Get slope of slopes:
    allslopes=np.array(allslopes)
    #    coeff,d=curve_fit(lin_fun,allslopes[:,0],allslopes[:,1])

#    ax.plot([0,0.7],[0,0.7],':k')
#        ax.plot([0,1.0],[coeff[1],coeff[0]*1.0+coeff[1]],':'+catcolors[icat],
#                label=cation+r': %1.2f$\frac{C_{A*}}{C_*}$%1.2f' %(coeff[0],coeff[1]))
    ax.plot([0,max(allslopes[:,0])],[0,max(allslopes[:,0])],'--k')
    ax.set_xlabel('C$_{*A}-C_{*}$ [e$^-$/V]')
    if seconderiv:
        ax.set_ylabel(r'$\frac{\partial^2 \Delta G[\phi]}{\partial \phi^2}$ [e$^-$/V]')
    else:
        ax.set_ylabel(r'$\left(\frac{(\partial \Delta G[\phi])}{\partial \phi}-1\right)\cdot C_*^{-1}$ [V]')
    fig.legend(loc='lower left',bbox_to_anchor=(0.75,0.2))
    fig.tight_layout()
    if seconderiv:
        fig.savefig(outdir+'dG_GC_second_deriv_v_capmismatch.pdf')
    else:
        fig.savefig(outdir+'dG_GC_vs_Capratio_all_cations.pdf')
    #plt.show()
    plt.close()


if __name__ == '__main__':
    main()

