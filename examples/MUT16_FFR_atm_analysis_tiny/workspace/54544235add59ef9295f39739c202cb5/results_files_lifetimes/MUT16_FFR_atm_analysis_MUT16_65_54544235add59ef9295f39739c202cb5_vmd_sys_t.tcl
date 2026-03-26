package require pbctools
color Display {Background} white
color scale method RWB
axes location off

# Load molecule and define selections
mol new ./MUT16_FFR_atm_analysis_MUT16_65_54544235add59ef9295f39739c202cb5_pi_clean.gro
mol addfile ./MUT16_FFR_atm_analysis_MUT16_65_54544235add59ef9295f39739c202cb5_full_pi_ref.xtc type {xtc} first 0 last -1 step 100 waitfor -1 
pbc box 

# Hide initial representation
mol delrep 0 top
set MUT16 [atomselect top "resid 1 to 172 or resid 173 to 344 or resid 345 to 516 or resid 517 to 688 or resid 689 to 860 or resid 861 to 1032 or resid 1033 to 1204 or resid 1205 to 1376 or resid 1377 to 1548 or resid 1549 to 1720 or resid 1721 to 1892 or resid 1893 to 2064 or resid 2065 to 2236 or resid 2237 to 2408 or resid 2409 to 2580 or resid 2581 to 2752 or resid 2753 to 2924 or resid 2925 to 3096 or resid 3097 to 3268 or resid 3269 to 3440 or resid 3441 to 3612 or resid 3613 to 3784 or resid 3785 to 3956 or resid 3957 to 4128 or resid 4129 to 4300 or resid 4301 to 4472 or resid 4473 to 4644 or resid 4645 to 4816 or resid 4817 to 4988 or resid 4989 to 5160 or resid 5161 to 5332 or resid 5333 to 5504 or resid 5505 to 5676 or resid 5677 to 5848 or resid 5849 to 6020 or resid 6021 to 6192 or resid 6193 to 6364 or resid 6365 to 6536 or resid 6537 to 6708 or resid 6709 to 6880 or resid 6881 to 7052 or resid 7053 to 7224 or resid 7225 to 7396 or resid 7397 to 7568 or resid 7569 to 7740 or resid 7741 to 7912 or resid 7913 to 8084 or resid 8085 to 8256 or resid 8257 to 8428 or resid 8429 to 8600 or resid 8601 to 8772 or resid 8773 to 8944 or resid 8945 to 9116 or resid 9117 to 9288 or resid 9289 to 9460 or resid 9461 to 9632 or resid 9633 to 9804 or resid 9805 to 9976 or resid 9977 to 10148 or resid 10149 to 10320 or resid 10321 to 10492 or resid 10493 to 10664 or resid 10665 to 10836 or resid 10837 to 11008 or resid 11009 to 11180"] ;

# Representation for MUT16
mol representation VDW 1.0 30.0
mol color ColorID 1
mol material Diffuse
mol selection "[$MUT16 text]"
mol addrep top
pbc box -color black
rotate y by 60
scale by 0.8

# --- Path to Binary (Adjust if necessary) ---
set tachyon_bin "/usr/local/vmd/lib/vmd/tachyon_LINUXAMD64"
set movie_file "MUT16_FFR_atm_analysis_MUT16_65_54544235add59ef9295f39739c202cb5_mol"
set num_frames [molinfo top get numframes]

# --- Lighting for Sharpness ---
display ambientocclusion on
display aoambient 0.5
display aodirect 0.8

# Loop through frames
for {set frame 0} {$frame < $num_frames} {incr frame} {
    animate goto $frame
    display update
    
    set dat_file "${movie_file}_${frame}.dat"
    set tga_file "${movie_file}_${frame}.tga"
    
    render Tachyon $dat_file
    
    if { [catch {exec $tachyon_bin -format TARGA -aasamples 12 -res 1920 1080 $dat_file -o $tga_file} msg] } {
        puts "Tachyon Error on frame $frame: $msg"
    }
}

# --- Movie Assembly (Lossless with Sharpening) ---
if { [file exists "${movie_file}_0.tga"] } {
    exec ffmpeg -framerate 24 -i "${movie_file}_%d.tga" \
        -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2,unsharp=5:5:1.0:5:5:0.0" \
        -c:v libx264 -crf 17 -preset slow -pix_fmt yuv420p -y "${movie_file}.mp4"
    
    exec ffmpeg -i "${movie_file}.mp4" -filter_complex "fps=15,scale=720:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -y "${movie_file}.gif"
    puts "RENDER SUCCESSFUL: ${movie_file}.mp4 generated."
}

# --- Final Cleanup ---
set files_to_delete [glob -nocomplain "${movie_file}*.tga" "${movie_file}*.dat"]
foreach file $files_to_delete {
    if {[file exists $file]} { file delete -force $file }
}
