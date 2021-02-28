#!/Strawberry/perl/bin/perl
use Getopt::Std;
use GD;
use GD::Text;
use GD::Text::Align;
require "./LipidomicsMapFunctions.pm";
use strict;

###### SYSTEM SETTINGS ######
my $FONT_PATH = "c:/Windows/Fonts"; # Ex. "c:/Windows/Fonts"

# Things still to do:
# I think the panels should be separately defined so we have the option of 
# putting them together in PS or gIMP to allow for more flexibility in layout
# This will work well for the fatty acid and sterols panels where there is a 
# mostly strictly defined structure. If there are too many FA types it could be
# problematic.
# For the sphingo and lipids panels, we can make the individual data panels and
# then small units that are the connectors, which could then be imported into
# PS or something for stitching together the final figure.

# We have a dataset that consists of a lipid compound string
# TYPE\(d?\d+:\d+(\/\d+:\d+)+(\(2OH\))?\)
# where the TYPE can be:
#   Cer
#   CL
#   DG
#   GalCer
#   PA
#   PC
#   PE
#   PG
#   PI
#   PI_Cer
#   PS
#   TG
# or some compounds without the parens
#   C\d+:\d+
#   Ergosta-7,22-dien-3-ol
#   Ergosterol
#   Campesterol
#   Ergost-7-en-3beta-ol
#   Lanosterol

# and an associated relative intensity value from 0.00 to 1.00

# The goal is to build an image for each type where the compounds are stacked by relative
# abundance, from highest to lowest, with a colored box, where the color indicates the
# relabund value with 1.0 being [red], 0.1 being [white] and 0.0 being [blue] (those colors
# will change to be color-blind friendly).

my $STYPE = {'Glycerolipids' => ['DG', 'TG'],
             'Glycerophospholipids' => ['PA', 'PC', 'PG', 'PS', 'LPG', 'CL', 'LPS', 'PI', 'PE', 'LPI', 'LPE', 'LPC'],
             'Sphingolipids' => ['PI_Cer', 'Cer(2OH)', 'HexCer', 'Cer'],
             'Sterols' => [],
             'Fatty acid' => ['FA']};

my $LIP = {}; # this will hold the lipid data, keyed by TYPE and CPD and holding the RELABUND
my $PROT = {}; # this will hold the proteomics data, keyed by STYPE
my %args; # this holds the command line arguments

# Handle input ################################################################
&getopts('l:p:m:o:f:h', \%args);
if ($args{'h'}) { &showHelp }
my $lipid_input = $args{'l'};
my $prot_input = $args{'p'};
my $output_file = $args{'o'};
my $output_format = $args{'f'};
my $FILL = {}; # this holds the fill color for each compound; relies on both the relative abundance data and the palettes

###############################################################################

# set up the output image #####################################################
my $xMargin = 16; # base x margin
my $yMargin = 16; # base y margin
my $pad = 3; # 3 px between elements
my $boxSize = 22; # size of elements
my $truecolor = 0;
# At some point, the image size could be determined by the input data, but for now we will have a preset size
my $imWidth = 2400;
my $imHeight = 2500;
my $im = new GD::Image($imWidth, $imHeight, $truecolor);
GD::Text->font_path($FONT_PATH) if ($FONT_PATH);
$im->setThickness(3); # set the base line thickness
# set up colors
our $palette = {'white' => $im->colorAllocate(255,255,255),
                'black' => $im->colorAllocate(0,0,0),
                'gray'  => $im->colorAllocate(127, 127, 127),
                'ltgray' => $im->colorAllocate(200, 200, 200),
                'Low' => $im->colorAllocate(251, 236, 94),
                'Moderate' => $im->colorAllocate(255, 168, 36),
                'High' => $im->colorAllocate(115, 74, 18)};
$palette->{'Very Low'} = $palette->{'white'};
$palette->{'Not Detected'} = $palette->{'gray'};
# This is a viridis palette
our $lipidPalette = &customPalette($im, {'R' =>  68, 'G' =>   2, 'B' => 84}, 
                                        {'R' =>  35, 'G' =>  138, 'B' => 141}, 
                                        {'R' => 251, 'G' => 231, 'B' => 37}
                                    );
#my $cividisPalette = &customPalette($im, {'R' => 0, 'G' => 32, 'B' => 77}, {'R'=> 124, 'G' => 124, 'B' => 120}, {'R' => 255, 'G' => 219, 'B' => 69})
#our $proteinPalette = {'Very Low' => $palette->{'white'}, 'Low' => $palette->{'Low'}, 'Moderate' => $palette->{'Moderate'}, 'High' => $palette->{'High'}};
###############################################################################

# draw the base map and fill fixed elements ###################################
# Read in the base map structure
my $MAP = &readMapStructure($args{'m'}, $palette, $im);
# some of the boxes in the base map are lipids or proteins, so read in that data
# to populate the colors for those boxes
# Read in the data to be displayed
$LIP = &readLipidData($lipid_input, $FILL);
$PROT = &readProteinData($prot_input, $FILL);

foreach my $layer (@{$MAP}) {
    foreach my $this (@{$layer}) {
        # this allows for overriding the default map box color with data from 
        # the lipid or protein data file.
        if ($this->{'shape'} eq 'e') {
            if (defined $FILL->{'e'}->{$this->{'label'}}) {
                $this->{'fill'} = $palette->{$FILL->{'e'}->{$this->{'label'}}};
            }   
        } elsif ($this->{'shape'} eq 'r') {
            if (defined $FILL->{'r'}->{$this->{'label'}}) {
                my $palidx = int($FILL->{'r'}->{$this->{'label'}} * 100);
                $this->{'fill'} = $lipidPalette->[$palidx];
            }
        }
        &drawObject($im, $this, $boxSize, $pad);
    }
}
###############################################################################

# fill in the keys ############################################################
&drawString($im, 
            {'color' => $palette->{'black'},
            'font' => "arial.ttf",
            'fontSize' => 24,
            'x' => 2135,
            'y' => 90,
            'label' => "Absolute protein\nabundance",
            'labJust' => 'center'});
my @protKey = ();
my $cx = 2065;
my $i = 0;
foreach my $key ('High', 'Moderate', 'Low', 'Very Low', 'Not Detected') {
    my $cy = 180 + $i++ * $boxSize;
    my $box_fill = $palette->{$key};
    $im->filledEllipse($cx, $cy, $boxSize, $boxSize, $box_fill);
    $im->ellipse($cx, $cy, $boxSize, $boxSize, $palette->{'black'});
    &drawString($im, {'x' => $cx + $boxSize + $pad, 
                      'y' => $cy - int($boxSize/2), 
                      'label' => $key, 
                      'labJust' => 'left',
                      'color' => $palette->{'black'}});
}

&drawString($im, 
            {'color' => $palette->{'black'},
            'font' => "arial.ttf",
            'fontSize' => 24,
            'x' => 2135,
            'y' => 362,
            'label' => "Relative\nlipid/metabolite\nabundance",
            'labJust' => 'center'});
my @lipidKey = (1, 0.5, 0.1, 0.05, 0, 'Not Detected');
my ($x1, $x2, $y1, $y2);
$x1 = 2065 - 0.5 * $boxSize;
$x2 = $x1 + $boxSize;
for (my $i=0; $i < 6; $i++) {
    my $y1 = $i * $boxSize + 465;
    my $y2 = $y1 + $boxSize;
    my $box_fill = $lipidKey[$i] eq "Not Detected" ? $palette->{'gray'} : $lipidPalette->[$lipidKey[$i] * 100]; 
    $im->filledRectangle($x1, $y1, $x2, $y2, $box_fill);
    $im->rectangle($x1, $y1, $x2, $y2, $palette->{'black'});
    &drawString($im, {'x' => $x2 + $boxSize/2 + $pad, 
                      'y' => $y1, 
                      'label' => $lipidKey[$i], 
                      'labJust' => 'left',
                      'color' => $palette->{"black"}});
}
###############################################################################

# define the container boxes for the compound lists (lipidomics data) #########
my %CONTAINER = ('Cer'    => {'shape' => 'r', 'x1' => 1824, 'y1' => 1692, 'x2' => 2114, 'y2' => 1820, 'fill' => 'white', 'color' => 'black'},
               'Cer(2OH)' => {'shape' => 'r', 'x1' => 1778, 'y1' => 1482, 'x2' => 2114, 'y2' => 1564, 'fill' => 'white', 'color' => 'black'},
               'CL'       => {'shape' => 'r', 'x1' => 1342, 'y1' => 1120, 'x2' => 1610, 'y2' => 1305, 'fill' => 'white', 'color' => 'black'},
               'DG'       => {'shape' => 'r', 'x1' =>  556, 'y1' =>  719, 'x2' =>  817, 'y2' => 1180, 'fill' => 'white', 'color' => 'black'},
               'FA'       => {'shape' => 'r', 'x1' => 1024, 'y1' =>  279, 'x2' => 1215, 'y2' =>  486, 'fill' => 'white', 'color' => 'black'},
               'HexCer'   => {'shape' => 'r', 'x1' => 1986, 'y1' => 1599, 'x2' => 2280, 'y2' => 1654, 'fill' => 'white', 'color' => 'black'},
               'LPC'      => {'shape' => 'r', 'x1' =>  171, 'y1' => 2152, 'x2' =>  368, 'y2' => 2360, 'fill' => 'white', 'color' => 'black'},
               'LPE'      => {'shape' => 'r', 'x1' => 1342, 'y1' => 2152, 'x2' => 1610, 'y2' => 2360, 'fill' => 'white', 'color' => 'black'},
               'LPG'      => {'shape' => 'r', 'x1' => 1342, 'y1' =>  975, 'x2' => 1610, 'y2' => 1082, 'fill' => 'white', 'color' => 'black'},
               'LPI'      => {'shape' => 'r', 'x1' => 1342, 'y1' => 2026, 'x2' => 1610, 'y2' => 2132, 'fill' => 'white', 'color' => 'black'},
               'LPS'      => {'shape' => 'r', 'x1' =>  949, 'y1' => 1373, 'x2' => 1164, 'y2' => 1480, 'fill' => 'white', 'color' => 'black'},
               'PA'       => {'shape' => 'r', 'x1' =>  949, 'y1' =>  719, 'x2' => 1164, 'y2' =>  903, 'fill' => 'white', 'color' => 'black'},
               'PC'       => {'shape' => 'r', 'x1' =>  554, 'y1' => 1266, 'x2' =>  797, 'y2' => 2360, 'fill' => 'white', 'color' => 'black'},
               'PE'       => {'shape' => 'r', 'x1' =>  949, 'y1' => 1494, 'x2' => 1164, 'y2' => 2360, 'fill' => 'white', 'color' => 'black'},
               'PG'       => {'shape' => 'r', 'x1' => 1342, 'y1' =>  719, 'x2' => 1610, 'y2' =>  903, 'fill' => 'white', 'color' => 'black'},
               'PI'       => {'shape' => 'r', 'x1' => 1342, 'y1' => 1404, 'x2' => 1610, 'y2' => 1910, 'fill' => 'white', 'color' => 'black'},
               'PI_Cer'   => {'shape' => 'r', 'x1' => 1824, 'y1' => 1279, 'x2' => 2114, 'y2' => 1437, 'fill' => 'white', 'color' => 'black'},
               'PS'       => {'shape' => 'r', 'x1' => 949,  'y1' => 1038, 'x2' => 1164, 'y2' => 1298, 'fill' => 'white', 'color' => 'black'},
               'TG'       => {'shape' => 'r', 'x1' => 171,  'y1' =>  719, 'x2' =>  416, 'y2' => 2092, 'fill' => 'white', 'color' => 'black'},
               );
###############################################################################

# iterate through the lipidomics data and draw the boxes ######################
while (my($stype, $types) = each %$STYPE) {
    foreach my $type (@$types) {
        #print "Making $type boxes...\n";
        # if (! defined $CONTAINER{$type}) { print "Pasta fa zool! Skipping!\n"; next; }
        my $this = $LIP->{$type};
        # draw a container box for each type
        &drawObject($im, $CONTAINER{$type});
        # make a side box for the type name
#        unless ($type eq "FA") {
        my $thisNAME = {%{$CONTAINER{$type}}};
        $thisNAME->{'x2'} = $thisNAME->{'x1'};
        $thisNAME->{'x1'} = $thisNAME->{'x1'} - 50;
        $thisNAME->{'fill'} = 'ltgray';
        &drawObject($im, $thisNAME);
        my $just = "center";
        my $label = $type eq "Cer(2OH)" ? "*" : $type eq "HexCer" ? "**" : $type;
        my (@bb) = &drawString($im, 
                                {'color' => $palette->{'black'},
                                'font' => "arial.ttf",
                                'fontSize' => 36,
                                'x' => $thisNAME->{'x1'} + $pad,
                                'y' => $thisNAME->{'y1'} + int(($thisNAME->{'y2'} - $thisNAME->{'y1'})/2),
                                'label' => $label,
                                'angle' => 90,
                                'labJust' => $just});
#       } # pairs with the unless ($type eq "FA")

        my $n = 0; # this is the count of elements to set the y-position
        foreach my $cpd (sort {$this->{$b}<=>$this->{$a}} keys %$this) {
            my $val = $this->{$cpd};
            
            # the x position of the element will be the margin
            my $x1 = $CONTAINER{$type}->{'x1'} + $xMargin;
            my $x2 = $x1 + $boxSize;
            # the y position will be the margin plus the (element size plus the padding) times the element count
            my $y1 = $CONTAINER{$type}->{'y1'} + $yMargin + ($n * ($boxSize + $pad));
            my $y2 = $y1 + $boxSize;
            # the fill color for the box will be based on the relabund value
            my $box_fill = $lipidPalette->[int($val * 100)];

            $im->filledRectangle($x1, $y1, $x2, $y2, $box_fill);
            $im->rectangle($x1, $y1, $x2, $y2, $palette->{'black'});

            # now to label the box
            #my @bounds = $im->stringFT($black, "arial.ttf", 14, 0, $x1 + $boxSize + $pad, $y1 + $boxSize, $cpd);
            my (@bb) = &drawString($im, 
                                   {'x' => $x1 + $boxSize + $pad,
                                    'y' => $y1,
                                    'label' => $cpd,
                                    'color' => $palette->{"black"}});
            $n++; # increment the count
            #print "$stype\t$type\t$cpd\t$val\n";
        }
    }
}

# write the image file
open(my $OUT, ">$output_file") or die "Can't open Lipid_map.png for write: $!\n";
binmode $OUT;
my $pngdata = $im->png();
print $OUT $pngdata;
close $OUT;

### Subroutines

sub readLipidData {
    my $file = shift;
    my $DATA = {};
    my $IN; # the input handle
    my $line; # the current line of the file being read

    if (-e $file && -r $file) {
        open($IN, $file) or die "Error reading $file: $!\n";
    } else {
        die "Cannot use $file; the file is either non-existant or unreadable.\n";
    }

    while ($line = <$IN>) {
        chomp $line;
        # my $type;
        my ($compound, $relabund, $type) = split(/\t/, $line, 3);
        if ($compound =~ /(\w+)\([td]?\d+:\d+(\/\d+:\d+)+(\(2OH\))?\)/) {
            $type = $1 if (! $type);
            if ($type =~ /P[CEGIS]/ && $compound =~ /\b0:0\b/) {
                $type = "L" . $type;
            }
            if ($type eq "Cer") { $type .= $3 } # this will append the (2OH) if it is present
            if ($type eq "GalCer") { $type = "HexCer" }
            $DATA->{$type}->{$compound} = $relabund;
        } elsif ($compound =~ /C\d+:\d+/) {
            $type = "FA";
            $DATA->{$type}->{$compound} = $relabund;
        } elsif ($compound =~ /sterol/i ||
                 $compound =~ /Ergost/i) {
            $type = "Sterol";
            $DATA->{$type}->{$compound} = $relabund;
            $FILL->{'r'}->{$compound} = $relabund;
        } else {
            $DATA->{$type}->{$compound} = $relabund;
            $FILL->{'r'}->{$compound} = $relabund;
#            print "What is this? $line\n";
        }
    }
    return($DATA);
}

sub readProteinData {
    my $file = shift;
    my $IN; # the input handle
    my $line; # the current line of the file being read

    if (-e $file && -r $file) {
        open($IN, $file) or die "Error reading $file: $!\n";
    } else {
        die "Cannot use $file; the file is either non-existant or unreadable.\n";
    }

    while ($line = <$IN>) {
        chomp $line;
        my ($gene, $abund, undef) = split(/\t/, $line, 3);
        $FILL->{'e'}->{$gene} = $abund;
    }
}

