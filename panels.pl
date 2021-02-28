#!/Strawberry/perl/bin/perl
use Getopt::Std;
use GD;
use GD::Text;
use GD::Text::Align;
use strict;
# might need to change to the full path for this to work
require "./LipidomicsMapFunctions.pm";

###### SYSTEM SETTINGS ######
my $FONT_PATH = ""; # Ex. "c:/Windows/Fonts"

# Handle input ################################################################
my %args; # this holds the command line arguments
getopts('hl:f:', \%args);
if ($args{'h'}) { &showHelp }
my $lipid_input = $args{'l'};
my $output_format = $args{'f'};

# color palette for the lipid data
our ($lipidPalette);

# defining the classes and types of lipids being tracked
my $STYPE = {'Glycerolipids' => ['DG', 'TG'],
             'Glycerophospholipids' => ['PA', 'PC', 'PG', 'PS', 'LPG', 'CL', 'LPS', 'PI', 'PE', 'LPI', 'LPE', 'LPC'],
             'Sphingolipids' => ['PI_Cer', 'Cer(2OH)', 'HexCer', 'Cer'],
             'Sterols' => [],
             'Fatty acid' => ['FA']
             };

# Read in the data to be displayed
my $LIP = &readLipidData($lipid_input);

my $xMargin = 5; # base x margin
my $yMargin = 5; # base y margin
my $pad = 3; # 3 px between elements
my $boxSize = 22;
my $truecolor = 0;
# this really should be a dynamic value that scales with the max label length for the type.
my $imWidth = 300;

while (my($class, $types) = each %$STYPE) {
    foreach my $type (@$types) {
        print "Making $type boxes...\n";
        my $this = $LIP->{$type};

        # calculate image height based on number of elements
        my $imHeight = ($boxSize + $pad) * scalar(keys %$this) + 2 * $yMargin;
        # instantiate the image
        my $im = new GD::Image($imWidth, $imHeight, $truecolor);
        GD::Text->font_path($FONT_PATH) if ($FONT_PATH);
        # define basic palette colors
        my $palette = {'white' => $im->colorAllocate(255,255,255),
                       'black' => $im->colorAllocate(0,0,0),
                       'gray'  => $im->colorAllocate(127, 127, 127),
                       'ltgray' => $im->colorAllocate(200, 200, 200),
                       'Low' => $im->colorAllocate(251, 236, 94),
                       'Moderate' => $im->colorAllocate(255, 168, 36),
                       'High' => $im->colorAllocate(115, 74, 18)};
        $palette->{'Very Low'} = $palette->{'white'};
        $palette->{'Not Detected'} = $palette->{'gray'};
        # this is replicating the cividis palette
        my $lipidPalette = &customPalette($im, {'R' =>  68, 'G' =>   2, 'B' => 84}, {'R' =>  35, 'G' =>  138, 'B' => 141}, {'R' => 251, 'G' => 231, 'B' => 37});

        # draw the panel
        $im->filledRectangle(51, 1, $imWidth, $imHeight, $palette->{"white"});
        $im->rectangle(50, 0, $imWidth - 1, $imHeight - 1, $palette->{"black"});
        # this is the panel label box
        $im->filledRectangle(1, 1, 51, $imHeight, $palette->{"ltgray"});
        $im->rectangle(0, 0, 50, $imHeight - 1, $palette->{"black"});
        my $just = "center";
        # the Ceramide panels tend to be too small for their labels, so we replace with asterisks
        my $label = $type eq "Cer(2OH)" ? "*" : $type eq "HexCer" ? "**" : $type;
        my ($x, $y) = (1 + $pad, 1 + int(($imHeight - 1)/2));
        my (@bb) = &drawString($im, 
                                {'color' => $palette->{"black"},
                                'fontSize' => 36,
                                'x' => $x,
                                'y' => $y,
                                'label' => $label,
                                'angle' => 90,
                                'labJust' => $just});

        # Now to draw the results boxes
        my $n = 0; # this is the count of elements to set the y-position
        foreach my $cpd (sort {$this->{$b}<=>$this->{$a}} keys %$this) {
            my $val = $this->{$cpd};
            
            # the x position of the element will be the margin plus the panel label width
            my $x1 = 50 + $xMargin;
            my $x2 = $x1 + $boxSize;
            # the y position will be the margin plus the (element size plus the padding) times the element count
            my $y1 = $yMargin + ($n * ($boxSize + $pad));
            my $y2 = $y1 + $boxSize;
            # the fill color for the box will be based on the relabund value
            my $box_fill = $lipidPalette->[int($val * 100)];
            $im->filledRectangle($x1, $y1, $x2, $y2, $box_fill);
            $im->rectangle($x1, $y1, $x2, $y2, 1);
            # now to label the box
            my (@bb) = &drawString($im, 
                                   {'x' => $x1 + $boxSize + $pad,
                                    'y' => $y1,
                                    'label' => $cpd},
                                    $palette->{"black"});
            $n++; # increment the count
        }

        # write the image file
        open(my $OUT, ">$type.png") or die "Can't open $type.png for write: $!\n";
        binmode $OUT;
        my $pngdata = $im->png();
        print $OUT $pngdata;
        close $OUT;
    }
}

