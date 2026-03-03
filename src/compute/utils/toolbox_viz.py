#create tcl files
def calc_vmd_labels(df,bool_dom,l_domains):
    str_or=' or '
    str_nl='\n'
    l_mol=[]
    l_beta=[]
    btb_grp=[]
    back_grp=[]
    idr1_grp=[]
    idr2_grp=[]
    val=40
    set_btb=''
    set_back=''
    set_idr1=''
    set_idr2=''
    selection_names = []
    #l_domains=['BTB', 'BACK']
    #d=0 #fummelfaktor
    for prot, df_prot in df.groupby('prot'):
        mol_grp=[]
        k=1/len(df.groupby('prot'))*20
        val+=k
        for index, row in df_prot.iterrows():
            mol_grp.append('resid {} to {}'.format(row['min'], row['max']))
            try:
                
                back_grp.append('resid {} to {}'.format(row['BACK_min'], row['BACK_max']))
            except:
                print('no such domain BACK')
                pass
            try:
                btb_grp.append('resid {} to {}'.format(row['BTB_min'], row['BTB_max']))
            except:
                print('no such domain BTB')
                pass
            try:
                idr1_grp.append('resid {} to {}'.format(row['IDR1_min'], row['IDR1_max']))
                
            except:
                print('no such domain IDR1')
                pass

            try:
                idr2_grp.append('resid {} to {}'.format(row['IDR2_min'], row['IDR2_max']))
            except:
                print('no such domain IDR2')
                pass
        set_mol='set {} [atomselect top "{}"] ;'.format(row['prot'],str_or.join(mol_grp))
        l_mol.append(set_mol)
        l_beta.append('${} set beta {}'.format(row['prot'],val))
        selection_names.append(row['prot'])
        #d=d+1
    if bool_dom:
        try:    
            set_btb='set {} [atomselect top "{}"] ;'.format('BTB',str_or.join(btb_grp))+str_nl
            set_back='set {} [atomselect top "{}"] ;'.format('BACK',str_or.join(back_grp))+str_nl
            l_beta.append('${} set beta {}'.format('BTB',100-k))
            l_beta.append('${} set beta {}'.format('BACK',k))
            
            selection_names.append('BTB')
            selection_names.append('BACK')

        except:
            print('no such domain BTB/BACK')
            pass
    
        try:    
            set_idr1='set {} [atomselect top "{}"] ;'.format('IDR1',str_or.join(idr1_grp))+str_nl
            

            if "resid" in set_idr1:
                selection_names.append('IDR1')
                l_beta.append('${} set beta {}'.format('IDR1',100-k))

            
        except:
            print('no such domain IDR1')

        try:
            set_idr2='set {} [atomselect top "{}"] ;'.format('IDR2',str_or.join(idr2_grp))+str_nl
            if "resid" in set_idr2:
                selection_names.append('IDR2')
                l_beta.append('${} set beta {}'.format('IDR2',k))
        except:
            print('no such domain IDR2')
            pass

            
    set_mols=str_nl.join(l_mol) +str_nl
    set_betas=str_nl.join(l_beta) +str_nl

    #add representations
    l_reps = []
    color_map = {
        'BTB':  4,
        'BACK': 16,
        'PEI1': 3,
        'PEI2': 22,
        'IDR1': 8,
        'IDR2': 8
        # Add other specific protein names and colors as needed
    }
    default_color = 30 # Fallback color for names not in the map

    # Get unique names while preserving order to avoid duplicate representations
    unique_names = list(dict.fromkeys(selection_names)) 
    for i, name in enumerate(unique_names):
        color = color_map.get(name, default_color)
        if name in l_domains:
            rep_block = f"""
    # Representation for {name}
    mol representation VDW
    mol color ColorID "{color}"
    mol material Transparent
    mol selection "[${name} text]"
    mol addrep top"""
        else:
             rep_block = f"""
    # Representation for {name}
    mol representation VDW
    mol color ColorID "{color}"
    mol material Diffuse
    mol selection "[${name} text]"
    mol addrep top"""
        l_reps.append(rep_block)
    set_reps = '\n'.join(l_reps)

    l_doms=[set_btb, set_back, set_idr1, set_idr2]
    return set_mols, l_doms,set_betas, set_reps


def create_vmd_file(domains_sys, top, traj='no', lab='no', type='NewCartoon 0.300000 10.000000 4.100000 0', step=1, nam="", bool_dom=True, l_domains=["BTB","BACK"]):
    #VMD Labeling
    #write domains to VMD file tcl
    #if bool_dom==True:
    mols, l_doms,betas,reps =calc_vmd_labels(domains_sys,bool_dom,l_domains)
    print("mols", mols,l_doms,betas,reps)
    #traj='mol addfile {{}} type {{xtc}} first 0 last -1 step 1 waitfor 1 0 '.format(traj)+'\n'
    mov = """
    # Define the output movie file name and format
    set movie_file "{}"
    set movie_format "tga" ;
    set output_format "mp4"
    
    # Set up the movie parameters
    set num_frames [molinfo top get numframes]
    set start_frame 0
    set end_frame [expr $num_frames - 1]
    
    # Loop through frames to generate .ps files
    for {{set frame $start_frame}} {{$frame <= $end_frame}} {{incr frame}} {{
        animate goto $frame
        set tga_file "${{movie_file}}_${{frame}}.ps"
        set extra_file "${{movie_file}}_${{frame}}.dat"
        render PostScript $tga_file
        render Tachyon $extra_file
    }}
    
    # Loop through .ps files to convert them to .tga format
    for {{set frame $start_frame}} {{$frame <= $end_frame}} {{incr frame}} {{
        set dat_file "${{movie_file}}_${{frame}}.ps"
        set tga_file "${{movie_file}}_${{frame}}.tga"
        # Convert .ps to .tga using an external Tachyon call
        convert $dat_file $tga_file
    }}
    
    # Combine frames into a movie
    exec ffmpeg -framerate 10 -i "${{movie_file}}_%d.tga" -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -pix_fmt yuv420p "${{movie_file}}.mp4"
    exec ffmpeg -i "${{movie_file}}.mp4" -c:v gif "${{movie_file}}.gif"
    
    # Cleanup intermediate frames
    set files_to_delete [glob -nocomplain "${{movie_file}}*.tga"]
    foreach file ${{files_to_delete}} {{
        if {{[file exists $file]}} {{
            file delete $file
        }}
    }}
    set files_to_delete2 [glob -nocomplain "${{movie_file}}*.ps"]
    foreach file ${{files_to_delete2}} {{
        if {{[file exists $file]}} {{
            file delete $file
        }}
    }}

    """.format(lab)

    content_end= """
    
    # update New representations
    #mol representation {}
    #mol color beta
    #mol color ColorID 0
    #mol color Chain
    #mol addrep top
    
    pbc box -color black
    rotate z by 0
    rotate x to 0
    rotate y by 60
    scale by 0.8
    """.format(type) 

    
    content = """
    package require pbctools
    color Display {{Background}} white
    color scale method RWB
    axes location off
    
    
    # Load molecule and define selections
    mol new {}  ;# path to your molecule file
    mol addfile {} type {{xtc}} first 0 last -1 step {} waitfor -1 
    pbc box 
    
    # Hide all representations
    mol delrep 0 top
    
    # Define domains
    """.format(top, traj,step)
    #if traj!='no':
    #    content +=traj

    content +=mols
    if bool_dom==True:
        for item in l_doms:
            print(item)
            if "resid" not in item:
                print('nono')
                continue
            content+=item
        #content +=btb
        #content +=back
        #content +=idr1
        #content +=idr2
        #content +=betas
    content += reps
    content +=content_end
    if lab!='no':
        content +=mov
    
    return content


def toggle_comment_lines(file_path, output_path, keywords, comment=True):
    """
    Comments or uncomments lines containing specific keywords in a file.

    :param file_path: Path to the input file.
    :param output_path: Path to the output file.
    :param keywords: List of keywords to look for in the lines.
    :param comment: If True, comments the lines; if False, uncomments them.
    """
    with open(file_path, 'r') as file, open(output_path, 'w') as output_file:
        for line in file:
            # Check if the line contains any of the specified keywords
            if any(keyword in line for keyword in keywords):
                if comment:
                    # If not already commented, comment the line
                    if not line.strip().startswith('#'):
                        line = '# ' + line
                else:
                    # If it's commented, uncomment the line
                    if line.strip().startswith('#'):
                        line = line.lstrip('#').strip() + '\n'
            # Write the line to the output file
            output_file.write(line)