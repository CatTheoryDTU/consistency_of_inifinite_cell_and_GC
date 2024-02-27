from scripts_analysis.general_parameters import *
from scripts_analysis.parse_tools import collect_results
from scripts_analysis.plot_tools import *
#from local_tools import    plot_charge_on_Na


def main():
    if 1:
        alldata,ceresults = collect_results(sizes)
    else:
        alldata,ceresults = pkl.load(open('parsed_data.pkl','rb'))

    create_dE_v_pot_plot(alldata)
#    plot_slope_vs_elneg(alldata)
#    plot_ce_dpot_plot(ceresults,alldata)
#    plot_GC_dne_vs_pot(alldata)
#    plot_dipoles_from_mean_potential(alldata)
#    plot_GC_capacitances(alldata=alldata)
#    plot_GC_dEcanondphi_vs_cap(alldata)
#    plot_dn_and_E_vs_Capmismatch(alldata)
#    plot_GC_dG_vs_Capratio(alldata)

#    plot_charge_on_Na() # Probably not working, quite old
main()
