def calc_vmd_labels(df, bool_dom, l_domains):
    str_or = ' or '
    str_nl = '\n'
    l_mol = []
    l_beta = []
    selection_names = []
    
    # Dictionary to hold lists of resid strings for each domain dynamically
    domain_groups = {dom: [] for dom in l_domains}
    
    val = 40
    num_prots = len(df.groupby('prot'))
    k = (1 / num_prots * 20) if num_prots > 0 else 0

    for prot, df_prot in df.groupby('prot'):
        mol_grp = []
        val += k
        
        for _, row in df_prot.iterrows():
            # 1. Handle Main Protein Selection
            mol_grp.append(f"resid {int(row['min'])} to {int(row['max'])}")
            
            # 2. Handle Arbitrary Domains if bool_dom is True
            if bool_dom:
                for dom in l_domains:
                    min_col, max_col = f"{dom}_min", f"{dom}_max"
                    if min_col in row and max_col in row:
                        if not (pd.isna(row[min_col]) or pd.isna(row[max_col])):
                            domain_groups[dom].append(f"resid {int(row[min_col])} to {int(row[max_col])}")
        
        prot_name = str(prot)
        l_mol.append(f'set {prot_name} [atomselect top "{str_or.join(mol_grp)}"] ;')
        l_beta.append(f'${prot_name} set beta {val}')
        selection_names.append(prot_name)

    l_dom_scripts = []
    if bool_dom:
        for i, dom in enumerate(l_domains):
            if domain_groups[dom]:
                selection_names.append(dom)
                script = f'set {dom} [atomselect top "{str_or.join(domain_groups[dom])}"] ;'
                l_dom_scripts.append(script + str_nl)
                beta_val = 100 - k if i % 2 == 0 else k
                l_beta.append(f'${dom} set beta {beta_val}')
            else:
                l_dom_scripts.append('') 

    vmd_colors = [1, 4, 7, 11, 12, 16, 22, 26, 29, 30, 31] 
    unique_names = list(dict.fromkeys(selection_names))
    dynamic_color_map = {name: vmd_colors[i % len(vmd_colors)] for i, name in enumerate(unique_names)}

    l_reps = []
    for name in unique_names:
        color = dynamic_color_map[name]
        material = "Transparent" if name in l_domains else "Diffuse"
        
        rep_block = f"""
# Representation for {name}
mol representation VDW 1.0 30.0
mol color ColorID {color}
mol material {material}
mol selection "[${name} text]"
mol addrep top"""
        l_reps.append(rep_block)

    return (
        str_nl.join(l_mol) + str_nl, 
        l_dom_scripts, 
        str_nl.join(l_beta) + str_nl, 
        '\n'.join(l_reps)
    )

def create_vmd_file(domains_sys, top, traj='no', lab='no', type='VDW', step=1, bool_dom=True, l_domains=["BTB","BACK"]):
    mols, l_doms, betas, reps = calc_vmd_labels(domains_sys, bool_dom, l_domains)
    
    # Updated Movie Logic for Tachyon and Sharper FFmpeg Output
    mov = f"""
# --- Path to Binary (Adjust if necessary) ---
set tachyon_bin "/usr/local/vmd/lib/vmd/tachyon_LINUXAMD64"
set movie_file "{lab}"
set num_frames [molinfo top get numframes]

# --- Lighting for Sharpness ---
display ambientocclusion on
display aoambient 0.5
display aodirect 0.8

# Loop through frames
for {{set frame 0}} {{$frame < $num_frames}} {{incr frame}} {{
    animate goto $frame
    display update
    
    set dat_file "${{movie_file}}_${{frame}}.dat"
    set tga_file "${{movie_file}}_${{frame}}.tga"
    
    render Tachyon $dat_file
    
    if {{ [catch {{exec $tachyon_bin -format TARGA -aasamples 12 -res 1920 1080 $dat_file -o $tga_file}} msg] }} {{
        puts "Tachyon Error on frame $frame: $msg"
    }}
}}

# --- Movie Assembly (Lossless with Sharpening) ---
if {{ [file exists "${{movie_file}}_0.tga"] }} {{
    exec ffmpeg -framerate 24 -i "${{movie_file}}_%d.tga" \\
        -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2,unsharp=5:5:1.0:5:5:0.0" \\
        -c:v libx264 -crf 17 -preset slow -pix_fmt yuv420p -y "${{movie_file}}.mp4"
    
    exec ffmpeg -i "${{movie_file}}.mp4" -filter_complex "fps=15,scale=720:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -y "${{movie_file}}.gif"
    puts "RENDER SUCCESSFUL: ${{movie_file}}.mp4 generated."
}}

# --- Final Cleanup ---
set files_to_delete [glob -nocomplain "${{movie_file}}*.tga" "${{movie_file}}*.dat"]
foreach file $files_to_delete {{
    if {{[file exists $file]}} {{ file delete -force $file }}
}}
"""

    content_header = f"""package require pbctools
color Display {{Background}} white
color scale method RWB
axes location off

# Load molecule and define selections
mol new {top}
mol addfile {traj} type {{xtc}} first 0 last -1 step {step} waitfor -1 
pbc box 

# Hide initial representation
mol delrep 0 top
"""

    content_end = f"""
pbc box -color black
rotate y by 60
scale by 0.8
"""

    # Combine parts
    content = content_header + mols
    if bool_dom:
        for item in l_doms:
            if "resid" in item:
                content += item
    
    content += reps
    content += content_end
    
    if lab != 'no':
        content += mov
    
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
