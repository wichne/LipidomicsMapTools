# lipidomics_map project

The purpose of this project is to ease visualization of lipidomic and proteomic data in a 
metabolic map format. This is a pilot project that is currently designed around only 
producing images for a particular study, and thus the lipidomics_map program is not very
flexible, and will likely not produce satisfactory visualizations for other studies. The
panels program was developed to easily produce elements that can be imported into PhotoShop
or GIMP to build customized maps. 

There are two programs
1. panels.pl
2. lipidomics_map.pl

and a file of common subroutines
LipidomicsMapFunctions.pm

The panels program takes as input a tab-delimited data file where the columns are
1. compound label
2. floating-point relative abundance value (range: 0 - 1)

Only lipids detected should be listed in the input file.
The program outputs an image file (PNG) for each lipid type. The image is a vertical stack
heatmap in descending abundance order. The colorscale is based onthe cividis palette.

The lipidomics_map program takes as input
1. a lipid abundance data file as described above
2. a protein abundance data file (described below)
3. a map structure file (described below)
4. the name for the output file

The protein abundance file is similar to the lipid abundance file, other than the
abundance is categorized as 'Very Low', 'Low', 'Moderate', or 'High'. Only detected
proteins should be listed in the input file.

The map structure file is a tab-delimited file that defines the pixel positions
of all elements of the map to be created. Elements include rectangles, ellipses,
lines, arcs, and arrows. See the example structure file for additional details.

REQUIREMENTS
Perl v5
GD module
GD::Text module
GD::Text::Align module
Math::Trig module
Getopt::Std module

INSTALLATION
Download the .pl and .pm files
Edit the '#!' lines in each file to point to your installation of Perl.
The LipidomicsMapFunctions.pm should be in the same directory as the .pl files.
You might need to edit the 'require' line to have an explicit path to this file.
Edit the $FONT_PATH line in each .pl file to point to your font library (should
be 'C:/Windows/Fonts' for Windows, can vary for Linux, but common locations are
'/usr/lib/share/fonts' and '/usr/share/fonts'). I am unfamiliar with Macs. The 
program should run without this setting, but the default fonts may be an unfortunate 
size/aesthetic.

EXECUTION

panels.pl -l lipid_input.tsv

lipidomics_map.pl -l lipid_input.tsv -p protein_input.tsv -o output_file.png -m map_structure_definition.tsv

