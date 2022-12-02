import matplotlib.pyplot as plt
# plt.style.use(['science','grid'])

k = 4
def define_sizes(col_width, k, ratio):
    return {
            'col_width' :3.5,
            'ratio' : 2, # Width/height ratio
            'width' : col_width*k,
            'height' : col_width*k/ratio,
    }

single_spcs = define_sizes(3.5, k, 2)
double_spcs = define_sizes(3.5, k, 2)
plt.rcParams.update({'font.size': 8*k, 'text.usetex':True})


def plot_plres(plres, t, option='ppl', axs=None):
    if axs is None:
        axs = plt.subplot()
    if option == 'ppl':
        axs.plot(plres.zhat, plres.ppl, label='$t=4.1$', linewidth=1.5)
    elif option == 'vpl':
        axs.plot(plres.z, plres.vpl, label='$t=4.1$', linewidth=1.5)
    axs.autoscale(enable = True, axis='both', tight=True)
    axs.set_xlabel(r'$\hat{z}$ - Active Power [kW]')
    axs.set_ylabel(r'PPL($\hat{z}$)')
    axs.legend(loc='lower right', fontsize='x-small',frameon=True, fancybox=True)
    return axs

def plot_ci(ymean, ci, axs=None, color=None):
    if axs is None:
        axs = plt.subplot()
    line, = axs.plot(ymean, color=color)
    axs.fill_between(ci.index, ci.zmax, ci.zmin, alpha=0.5,color=line.get_color(), linewidth=0)
    axs.autoscale(enable = True, axis='both', tight=True)
    axs.set_xlabel('Time (s)')
    axs.set_ylabel('Active Power (kW)')
    return axs


def plot_ci_size(ci, axs=None):
    if axs is None:
        axs = plt.subplot()
    axs.plot(ci.max - ci.min)
    axs.autoscale(enable = True, axis='both', tight=True)
    axs.set_xlabel('Time (s)')
    return axs