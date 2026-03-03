import os
import re
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sbn


def flatten_mi(labels):
    """Turn MultiIndex tuples into one-line strings like 'A | 12 | foo'."""
    return [" | ".join(map(str, lab)) if isinstance(lab, tuple) else str(lab)
            for lab in labels]

def decimate(n, target=60):
    """Show at most `target` ticks by computing a stride."""
    import math
    return max(1, math.ceil(n / target))




def sliding_window_mean(data, window_size):

    num_windows = len(data) // window_size
    smoothed_data = np.zeros(num_windows)
    smoothed_index = np.zeros(num_windows, dtype=int)
    
    for i in range(num_windows):
        start_index = i * window_size
        end_index = min((i + 1) * window_size, len(data))
        window_data = data[start_index:end_index]
        smoothed_data[i] = np.mean(window_data)
        smoothed_index[i] = start_index
    
    return smoothed_data, smoothed_index

def load_data(file_path):
    x = []
    y = []

    with open(file_path, "r") as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith("#") or line.startswith("@"):
                continue
            values = line.split()
            x.append(float(values[0]))
            y.append(float(values[1]))
            print("loaded data: ",len(x), len(y))

    return {'x': x, 'y': y}

def plot_data(x_in, y_in, x_label, y_label, title, plots_folder, axis_label_size=9,avrg=False, window=2, t_step=1, every_x=1, rot=0, lim=False,llim_x=0, ulim_x=0,llim_y=0,ulim_y=0):
    # Extract the x and y values from the data
    if avrg==True:
            mean_dat, mean_ind = sliding_window_mean(y_in, window)
            #print("means", len(mean_dat), len(mean_ind))
    # Generate the plot
    plt.figure(figsize=(8, 6))
    plt.plot(x_in, y_in, '+-')
    if avrg==True:
        #plt.plot(mean_ind*t_step*np.ones(len(mean_ind)), mean_dat,'+-', label='window mean'+str(window)+' win')
        plt.plot([x * t_step for x in mean_ind], mean_dat,'+-', label='window mean'+str(window)+' win '+ str(t_step) + 'factor') 

    if lim=='x':
        plt.xlim(llim_x,ulim_x)
        print("set x lim")
    elif lim=='y':
        plt.ylim(llim_y,ulim_y)
        print("set y lim")
        #print("Joint plot y-limits:", joint_ax.ylim())
    elif lim=='xy':
        plt.xlim(llim_x,ulim_x)
        plt.ylim(llim_y,ulim_y)
        print("set x lim to ",llim_x,ulim_x)
        print("set y lim to ",llim_y,ulim_y)
    
    plt.xlabel(x_label, fontsize=axis_label_size)
    plt.ylabel(y_label, fontsize=axis_label_size)
    min_tot=min(x_in)
    max_tot=max(x_in)
    plt.xticks(np.arange(min_tot, max_tot + 1, int(t_step)*every_x),rotation=rot)
    plt.title(title, fontsize=axis_label_size + 2)

    # Regular expression pattern to match spaces, slashes, brackets, and periods
    pattern = r"[ /()\[\]{}.]"

    # Replace matched characters with an underscore
    title_string = re.sub(pattern, '_', title)
    
    plot_filename = f'{title_string}.png'
    plot_filepath = os.path.join(plots_folder, plot_filename)
    plt.savefig(plot_filepath)


def generate_and_save_plots(folder_path, data_file_names, x_labels, y_labels, titles, y_joint, joint_plot='x', axis_label_size=9, concate=False, lim=False, llim_x=0, ulim_x=0,llim_y=0,ulim_y=0, mean=False, win=1, step_to_time=[1], enumerate_x=False, every_x=1, rot=0, save='folder_path', dat=False, n_ticks=30):

    
    # Create 'plots' subfolder if it doesn't exist
    if save=='folder_path':
        plots_folder = os.path.join(folder_path, 'plots')
    else:
        plots_folder=save
        plots_folder = os.path.join(plots_folder, 'plots')
        
    os.makedirs(plots_folder, exist_ok=True)
    #print('check plots folder', plots_folder)

    if not joint_plot:
        for idx, (data_file_name, x_label, y_label, title, step_to_time) in enumerate(zip(data_file_names, x_labels, y_labels, titles, step_to_time)):

            if dat==True:
                #print(idx)
                data=dict({"k":0})
                data['x']=list(data_file_name[0])#{'x':data_file_name[0], 'y':data_file_name[1]}
                data['y']=list(data_file_name[1])
            else:
                data = load_data(os.path.join(folder_path, data_file_name))
                #print('plotted data',len(data['x']), len(data['y']))

            plot_data(data['x'], data['y'], x_label, y_label, title, plots_folder, avrg=mean, window=win, t_step=step_to_time, lim=lim, llim_x=llim_x, ulim_x=ulim_x,llim_y=llim_y,ulim_y=ulim_y)
        
    if joint_plot:
        x_max=0  
        min_tot=0
        max_tot=0
        y_values_all = []
        # Generate joint plots and save them to a single file
        joint_fig, joint_ax = plt.subplots()
        #for idx, (data_file_name, x_label, y_label, title) in enumerate(zip(data_file_names, x_labels, y_labels, titles)):
        for idx, (data_file_name,x_label, y_label, title) in enumerate(zip(data_file_names, x_labels, y_labels, titles)):

            if dat==True:
                #print(idx)
                data=dict({"k":0})
                data['x']=list(data_file_name[0])#{'x':data_file_name[0], 'y':data_file_name[1]}
                data['y']=list(data_file_name[1])
            else:
                data = load_data(os.path.join(folder_path, data_file_name))
                
            x = data['x']
            
            if enumerate_x:
                x = [i+1 for i in np.arange(len(data['x']))] 
            if concate==True:
                #x = data['x']+x_max*np.ones(len(data['x']))
                x = x+x_max*np.ones(len(x))
                #print("x_axis_values", x)
            y = data['y']
            #joint_ax.scatter(x, y, label=title)
            y_values_all.append(y)
            joint_ax.plot(x, y, '+-', label=title)
            
            x_max=x_max+max(data['x'])
            min_tot=min(min_tot,min(x))
            max_tot=max(max_tot,max(x))
        try:    
            mean_y = np.mean(np.array(y_values_all), axis=0)    
            joint_ax.plot(x, mean_y, 'b--', label='Mean', linewidth=2)   
        except:
            pass
        
        joint_ax.set_xlabel(x_labels[0], fontsize=axis_label_size)
        tick_step = max(1, int((max_tot - min_tot) / n_ticks))  # Ensure at least step of 1
        plt.xticks(np.arange(min_tot, max_tot + 1, step=tick_step), rotation=rot)

        
        print("max_step_to_time", max(step_to_time))
        #does not work with only one axis set and the other one automatically scaled
        if lim=='x':
            joint_ax.set_xlim(llim_x,ulim_x)
            print("set x lim")
        elif lim=='y':
            joint_ax.set_ylim(llim_y,ulim_y)
            print("set y lim")
            print("Joint plot y-limits:", joint_ax.get_ylim())
        elif lim=='xy':
            joint_ax.set_xlim(llim_x,ulim_x)
            joint_ax.set_ylim(llim_y,ulim_y)
            print("set x lim to ",llim_x,ulim_x)
            print("set y lim to ",llim_y,ulim_y)
        joint_ax.set_ylabel(y_joint, fontsize=axis_label_size)
        joint_ax.set_title(y_joint+ ' overview', fontsize=axis_label_size + 2)
        joint_ax.tick_params(labelsize=axis_label_size)
        
        #joint_ax.legend()
        joint_ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), borderaxespad=0.)
        joint_fig.savefig(os.path.join(plots_folder, y_joint+"_joint_plots.png"),bbox_inches='tight')


def reorder_dataframe(dat, o, attr):
    """Reorder MultiIndex dataframe based on categorical order tuples (o)."""
    if attr != "res_org":
        a1, a2 = o[0], o[1]
        pairs = [(x, y) for x in a1 for y in a2]
    else:
        a1 = o[0]
        a2 = dat.index.get_level_values(1)
        pairs = [(x, y) for x in a1 for y in a2]

    ind_filt = [p for p in pairs if p in dat.index]
    col_filt = [p for p in pairs if p in dat.columns]

    dat = dat.reindex(index=ind_filt)
    dat = dat.reindex(columns=col_filt)
    dat = dat.fillna(0.0).astype("float64")
    return dat


def validate_data(dat):
    """Raise errors for empty or invalid dataframes."""
    if dat.size == 0:
        raise ValueError("Empty after ordering; check MultiIndex structure.")
    if np.isnan(dat.to_numpy()).all():
        raise ValueError("All-NaN after conversion; nothing to plot.")


def plot_heatmap(
    dat,title,filename, target_ticks=60, figsize=(40, 40), dpi=300, show_colorbar=True,ftnsz=10
,values=False,min=0):
    """Render and save a single heatmap with decimated tick labels."""
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    ax.set_aspect("equal")

    hm = sbn.heatmap(
        dat,
        xticklabels=False,
        yticklabels=False,
        square=False,
        cbar=show_colorbar,
        cbar_kws={"shrink": 0.5} if show_colorbar else None,
        cmap="Spectral_r",
        annot=values,
        annot_kws={"size": 1.5*ftnsz},
        fmt=".6f",
        #linewidths=0.1,
        #linecolor="black",
        vmin=min,
        vmax=1,
        ax=ax
    )
    hm.invert_yaxis()

    #if show_colorbar:
    #    hm.collections[0].colorbar.set_label(
    #        "p(i,j)|p(l,m) - conditional probability", fontsize=2*ftnsz
    #        
    #    )
    #    cbar.ax.tick_params(labelsize=1.8 * ftnsz) 
    if show_colorbar:
        cbar = hm.collections[0].colorbar
        cbar.set_label(
            "p(i,j)|p(l,m) - conditional probability",
            fontsize=2 * ftnsz,
            labelpad=10,
        )
        cbar.ax.tick_params(labelsize=1.8 * ftnsz) 


    
    # Tick decimation
    x_labels = flatten_mi(dat.columns)
    y_labels = flatten_mi(dat.index)
    step_x = decimate(len(x_labels), target=target_ticks)
    step_y = decimate(len(y_labels), target=target_ticks)
    x_idx = np.arange(0, len(x_labels), step_x)
    y_idx = np.arange(0, len(y_labels), step_y)
    ax.set_xticks(x_idx + 0.5)
    ax.set_yticks(y_idx + 0.5)
    ax.set_xticklabels(
        [x_labels[i] for i in x_idx], rotation=90, ha="center", fontsize=2*ftnsz
    )
    ax.set_yticklabels(
        [y_labels[i] for i in y_idx], rotation=0, ha="right", fontsize=2*ftnsz
    )
    ax.tick_params(axis="x", pad=40)
    ax.tick_params(axis="y", pad=40)
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.title(title, fontsize=2*ftnsz)
    fig.savefig(filename, dpi=dpi)
    plt.close(fig)


def plot_subblocks(dat, label, r0, c0, attr, query_n, results_path, file_res=300,ftnsz=10,values=False):
    """Plot each first-level block separately and also a large overview grid."""
    row_groups = dat.index.get_level_values(0).unique()
    col_groups = dat.columns.get_level_values(0).unique()

    n_rows, n_cols = len(row_groups), len(col_groups)
    fig_big, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(n_cols * 10, n_rows * 10),
        constrained_layout=False,
        squeeze=False,
    )

    for i, rg in enumerate(row_groups):
        for j, cg in enumerate(col_groups):
            ax = axes[i, j]
            sub_dat = dat.loc[
                dat.index.get_level_values(0) == rg,
                dat.columns.get_level_values(0) == cg,
            ]
            if sub_dat.empty:
                ax.axis("off")
                continue



            sbn.heatmap(
                sub_dat,
                xticklabels=True,
                yticklabels=True,
                square=False,
                cbar=False,
                cmap="Spectral_r",
                annot=values,
                #linewidths=0.05,
                #linecolor="black",
                fmt=".4f",
                annot_kws={"size": 2*ftnsz},
                vmin=0,
                vmax=1,
                ax=ax,
            ).invert_yaxis()
            ax.set_title(f"{rg} - {cg}", fontsize=ftnsz)
            ax.set_aspect("equal")
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.tick_params(axis="x", labelsize=ftnsz,pad=12)
            ax.tick_params(axis="y", labelsize=ftnsz,pad=12)
            x_labels = [str(x) for x in sub_dat.columns.get_level_values(1)]
            y_labels = [str(y) for y in sub_dat.index.get_level_values(1)]
            ax.set_xticks(np.arange(len(x_labels)) + 0.5)
            ax.set_yticks(np.arange(len(y_labels)) + 0.5)
            ax.set_xticklabels(x_labels, rotation=90)#, ha="center", fontsize=20)
            ax.set_yticklabels(y_labels, rotation=0)#, ha="right", fontsize=20)

            # Individual file for each block
            subfile = f"{results_path}/{label}_given_{r0}_{c0}_plot_{attr}_by_{query_n}_{rg}_to_{cg}.png"
            print(subfile)
            plot_heatmap(
                sub_dat,
                title=f"{r0}_{c0}_{attr}_by_{query_n} | {rg}-{cg}",
                filename=subfile,
                target_ticks=30,
                figsize=(25, 25),
                dpi=file_res,
            )
    cbar_ax = fig_big.add_axes([0.92, 0.25, 0.012, 0.4])  # narrower + shorter
    norm = plt.Normalize(vmin=0, vmax=1)
    sm = plt.cm.ScalarMappable(cmap="Spectral_r", norm=norm)
    sm.set_array([])
    
    cbar = fig_big.colorbar(
        sm,
        cax=cbar_ax,
        label="p(i,j)|p(l,m) - conditional probability"
    )
    
    # 🔧 Scale down the colorbar aesthetics
    cbar.ax.tick_params(labelsize=ftnsz)       # smaller tick labels
    cbar.set_label(
        "p(i,j)|p(l,m) - conditional probability",
        fontsize=ftnsz,
        labelpad=ftnsz
    )
    fig_big.subplots_adjust(
    left=0.05,
    right=0.87,   # leave room for colorbar
    bottom=0.05,
    top=0.95,
    wspace=0.2,   # ⬅ increase horizontal spacing
    hspace=0.2    # ⬅ increase vertical spacing
)

    
    # Shared colorbar for the overview
    #cbar_ax = fig_big.add_axes([0.93, 0.15, 0.02, 0.7])
    #norm = plt.Normalize(vmin=0, vmax=1)
    #sm = plt.cm.ScalarMappable(cmap="Spectral_r", norm=norm)
    #sm.set_array([])
    #fig_big.colorbar(sm, cax=cbar_ax, label="p(i,j)|p(l,m) - conditional probability")

    fig_big.suptitle(
        f"{r0}_{c0}_{attr}_by_{query_n})",
        fontsize=ftnsz,
    )
    overview_file = f"{results_path}/{label}_given_{r0}_{c0}_plot_{attr}_by_{query_n}_overview.png"
    print(overview_file)
    fig_big.savefig(overview_file, dpi=300)
    plt.close(fig_big)