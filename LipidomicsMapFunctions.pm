#!/Strawberry/perl/bin/perl

use Math::Trig;
use GD;
use GD::Text;
use GD::Text::Align;
use strict;

sub drawString {
    my $im = shift;
    my $obj = shift;
    my $color = shift;

    my %config = (font => 'arial',
                  ptsize => 14,
                  valign => 'top',
                  halign => 'left',
                  color => $color);
    my $angle = 0;
    

    if (defined $obj->{'valign'}) { $config{'valign'} = $obj->{'valign'}}
    if (defined $obj->{'labJust'}) {
        #print "Setting halign to " . $obj->{'labJust'} . " for " . $obj->{'label'} . "\n";
        $config{halign} = $obj->{'labJust'};
    }
    if (defined($obj->{'font'})) { $config{font} = $obj->{'font'} }
    if (defined($obj->{'fontSize'})) { $config{ptsize} = $obj->{'fontSize'} }
    if (defined($obj->{'color'})) { $config{color} = $obj->{'color'} }
    $config{text} = $obj->{'label'};

    if (defined($obj->{'angle'})) { $angle = rad($obj->{'angle'}) }
 
    #print STDOUT "\t" . $config{'text'} . "\n";
    my $aln = GD::Text::Align->new($im, %config);
    my @bb = $aln->draw($obj->{'x'}, $obj->{'y'}, $angle);
    # print join(" ", @bb) . "\n";
    # $im->rectangle($bb[0], $bb[1], $bb[4], $bb[5], $color);
}

sub customPalette {
    my $im = shift;
    my $c1 = shift;
    my $c2 = shift;
    my $c3 = shift;

    # To reserve space for other colors, we will only do a 101 color scale.
    # Also, since the scale gradation is 1 to 0.1 to 0, the one gradient will be 0 to 10 and the other will be 11 to 100

    my @palette;
    $palette[0]   = $im->colorAllocate($c1->{R}, $c1->{G}, $c1->{B});
    my $offset = 0;
    my $steps = 10;
    for (my $i=0; $i<$steps; $i++) {
        my $thisR = $c1->{R} - int(((($c1->{R} - $c2->{R})/$steps) * $i) + 0.5);
        my $thisG = $c1->{G} - int(((($c1->{G} - $c2->{G})/$steps) * $i) + 0.5);
        my $thisB = $c1->{B} - int(((($c1->{B} - $c2->{B})/$steps) * $i) + 0.5);
        $palette[$offset + $i] = $im->colorAllocate($thisR, $thisG, $thisB);
        # print "$i: $thisR, $thisG, $thisB\n";
    }
    $palette[$steps] = $im->colorAllocate($c2->{R}, $c2->{G}, $c2->{B});

    $offset = $steps;
    $steps = 91;
    for (my $i=0; $i<$steps; $i++) {
        my $thisR = $c2->{R} - int(((($c2->{R} - $c3->{R})/$steps) * $i) + 0.5);
        my $thisG = $c2->{G} - int(((($c2->{G} - $c3->{G})/$steps) * $i) + 0.5);
        my $thisB = $c2->{B} - int(((($c2->{B} - $c3->{B})/$steps) * $i) + 0.5);
        $palette[$offset+$i] = $im->colorAllocate($thisR, $thisG, $thisB);
        # print $i + $offset, ": $thisR, $thisG, $thisB\n";
    }
    return \@palette;
}
sub readLipidData {
    my $file = shift;
    my $FILL = shift;

    my $DATA = {};
    my $IN; # the input handle
    my $line; # the current line of the file being read
    if (-e $file && -r $file) {
        open($IN, $file) or die "Error reading $file: $!\n";
    } else {
        die "Cannot use lipid file $file; the file is either non-existant or unreadable.\n";
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
        } else {
            $DATA->{$type}->{$compound} = $relabund;
#            print "What is this? $line\n";
        }
        $FILL->{'r'}->{$compound} = $relabund;
    }
    return($DATA);
}

sub readProteinData {
    my $file = shift;
    my $FILL = shift; # this is the data structure that holds the relative abundance values
    my $IN; # the input handle
    my $line; # the current line of the file being read

    if (-e $file && -r $file) {
        open($IN, $file) or die "Error reading $file: $!\n";
    } else {
        die "Cannot use protein file $file; the file is either non-existant or unreadable.\n";
    }

    while ($line = <$IN>) {
        chomp $line;
        my ($gene, $abund, undef) = split(/\t/, $line, 3);
        $FILL->{'e'}->{$gene} = $abund;
    }
}

sub rad {
    my $deg = shift;
    return($deg * pi/180);
}

sub drawObject {
    my $im = shift;
    my $obj = shift;
    my $boxSize = shift;
    my $pad = shift;

    if (! $obj->{color}) { $obj->{color} = 1; }

    # now figure out what we're drawing and draw it.
    if ($obj->{'shape'} eq 'l') {
        $im->line($obj->{'x1'}, $obj->{'y1'}, $obj->{'x2'}, $obj->{'y2'}, $obj->{color});
    } elsif ($obj->{'shape'} eq 'e') {
#        print "drawing '" . $obj->{label} . "'\n";
        $im->filledEllipse($obj->{'cx'}, $obj->{'cy'}, $obj->{'width'}, $obj->{'height'}, $obj->{'fill'});
        $im->ellipse($obj->{'cx'}, $obj->{'cy'}, $obj->{'width'}, $obj->{'height'}, $obj->{color});
        my $eoffset = 0.5 * $boxSize;
        if ($obj->{'labPos'}) {
            my $valign = "top";
            my ($x, $y) = ($obj->{'cx'} - $eoffset, $obj->{'cy'} - $eoffset + $pad);
            if ($obj->{'labPos'} =~ /^([LRTB])([LRTB]?)$/) {
                if ($1 eq "R") { $x = $x + $boxSize + 2*$pad; }
                elsif ($1 eq "L") { $x = $x - $pad; }
                elsif ($1 eq "T") { $y = $y - $pad; $valign = 'bottom'; $x = $obj->{'cx'};}
                elsif ($1 eq "B") { $y = $y + $boxSize + 2*$pad; $x = $obj->{'cx'}; }

                if ($2 eq "R") { $x = $obj->{'cx'} + $eoffset + $pad; }
                elsif ($2 eq "L") { $x = $obj->{'cx'} - $eoffset - $pad; }
                elsif ($2 eq "T") { $y = $y - $eoffset; }
                elsif ($2 eq "B") { $y = $y + $eoffset; }
                
            } else { print "How to interpret " . $obj->{'labPos'} . "?\n"; }

            #my @bounds = $im->stringFT($black, "arial.ttf", 14, 0, $x, $y + $boxSize, $obj->{'label'});
            my (@bb) = &drawString($im, 
                        {'x' => $x,
                        'y' => $y,
                        'label' => $obj->{'label'},
                        'fontSize' => 14,
                        'valign' => $valign,
                        'labJust' => $obj->{'labJust'},
                        'color' => $obj->{'color'}
                        });

#            $im->string(gdLargeFont, $x, $y, $obj->{'label'}, $black);
        }
    } elsif ($obj->{'shape'} eq 'r') {
#        print "drawing " . $obj->{label} . "\n";
        $im->filledRectangle($obj->{'x1'}, $obj->{'y1'}, $obj->{'x2'}, $obj->{'y2'}, $obj->{fill});
        $im->rectangle($obj->{'x1'}, $obj->{'y1'}, $obj->{'x2'}, $obj->{'y2'}, $obj->{color});
        if ($obj->{'labPos'}) {
            my ($x, $y) = ($obj->{'x1'}, $obj->{'y1'} + $pad);
            my $valign = "top";
            if ($obj->{'labPos'} =~ /^([LRTB])([LRTB]?)$/) {
                if ($1 eq "R") { $x = $x + $boxSize + 2*$pad; }
                elsif ($1 eq "L") { $x = $x - $pad; }
                elsif ($1 eq "T") { $y = $y - $pad; $valign = 'bottom'; $x = $x + int(($obj->{'x2'} - $obj->{'x1'})/2);}
                elsif ($1 eq "B") { $y = $y + $boxSize + 2*$pad; $x = $x + int(($obj->{'x2'} - $obj->{'x1'})/2);}

                if ($2 eq "R") { $x = $obj->{'x2'} + $pad; }
                elsif ($2 eq "L") { $x = $obj->{'x1'} - $pad; }
                elsif ($2 eq "T") { $y = $y - 0.5 * $boxSize; }
                elsif ($2 eq "B") { $y = $y + 0.5 * $boxSize; }
                
            } else { print "How to interpret " . $obj->{'labPos'} . "?\n"; }

            #my @bounds = $im->stringFT($black, "arial.ttf", 14, 0, $x, $y + $boxSize, $obj->{'label'});
            my (@bb) = &drawString($im, 
                        {'x' => $x,
                         'y' => $y,
                         'label' => $obj->{'label'},
                         'fontSize' => 14,
                         'valign' => $valign,
                         'labJust' => $obj->{'labJust'},
                         'color' => $obj->{'color'}
                        });
        }
    } elsif ($obj->{'shape'} eq 'a') {
        $im->arc($obj->{'cx'}, $obj->{'cy'}, $obj->{'width'}, $obj->{'height'}, $obj->{'start'}, $obj->{'end'}, $obj->{color});
    } elsif ($obj->{'shape'} eq 't') {
        #print "Making an arrow!\n";
        &Arrow($im, $obj->{'ptx'}, $obj->{'pty'}, 10, 5, $obj->{'dir'}, $obj->{'color'});
    }
}

sub hexToNum {
    # this take a 6-character hex string and turns it into 3 numerics, you know, for colors
    my $hex = shift;
    my @num;
    if ($hex =~ /^([0-9a-f]{2})([0-9a-f]{2})([0-9a-f]{2})$/) {
        #print "$1, $2, $3 -> ";
        @num = (hex($1), hex($2), hex($3));
        #print join(", ", @num) . "\n";
    }
    return @num;
}

sub readMapStructure {
    my $file = shift;
    my $palette = shift; # hash keyed by color name storing palette index
    my $im = shift; # in case we need to define new colors
    my $MAP = [];
    open(my $IN, $file) or die "Can't open map file $file for reading: $!\n";
    while (my $line = <$IN>) {
        chomp($line);
        my ($depth,
            $label, 
            $shape, 
            $x1, 
            $y1, 
            $x2, 
            $y2, 
            $width, 
            $height, 
            $labPos, 
            $labJust, 
            $arrow, 
            $start, 
            $end,
            $fill,
            $color) = split(/\t/, $line);
        $label =~ s/[\"]//g;
#        if (defined $FILL->{$shape}->{$label}) {
#            print "Setting $label fill from FILL\n";
#            $fill = $FILL->{$shape}->{$label};
#        }
        if (! $fill) { 
            $fill = "gray";
        }
        if (! defined $palette->{$fill}) {
            my ($r, $g, $b) = &hexToNum($fill);
            #print "Ooh! got " . $obj->{'fill'} . " to be $r, $g, $b\n";
            $palette->{$fill} = $im->colorAllocate($r, $g, $b);
        }

        if (! $color) { 
            $color = "black";
        }
        if (! defined $palette->{$color}) {
            my ($r, $g, $b) = &hexToNum($color);
            #print "Ooh! got " . $obj->{'fill'} . " to be $r, $g, $b\n";
            $palette->{$color} = $im->colorAllocate($r, $g, $b);
        }

        my $this;
        if ($labJust) {
            if ($labJust =~ /^r/i) { $labJust = "right"}
            elsif ($labJust =~ /^l/i) { $labJust = "left"}
            elsif ($labJust =~ /^c/i) { $labJust = "center"}
            else { warn "What is up with $line:\nspecifically '$labJust'. Using 'left'.\n"; $labJust = "left";}
        }
        if ($shape eq "r") {
            $this = {'shape' => 'r',
                     'x1' => $x1,
                     'x2' => $x2,
                     'y1' => $y1,
                     'y2' => $y2,
                     'label' => $label,
                     'labPos' => $labPos,
                     'labJust' => $labJust};
            $this->{'fill'} = $palette->{$fill} if $fill;
            $this->{'color'} = $palette->{$color} if $color;
        } elsif ($shape eq "e") {
            if (! $height) { $height = $width};
            $this = {'shape' => 'e',
                        'cx' => $x1,
                        'cy' => $y1,
                        'width' => $width,
                        'height' => $height,
                        'label' => $label,
                        'labPos' => $labPos,
                        'labJust' => $labJust};
            $this->{'fill'} = $palette->{$fill} if $fill;
            $this->{'color'} = $palette->{$color} if $color;
        } elsif ($shape eq "l") {
            $this = {'shape' => 'l',
                        'x1' => $x1,
                        'x2' => $x2,
                        'y1' => $y1,
                        'y2' => $y2,
                        'arrow' => $arrow};
        } elsif ($shape eq "a") {
            $this = {'shape' => 'a',
                        'cx' => $x1,
                        'cy' => $y1,
                        'width' => $width,
                        'height' => $height,
                        'start' => $start,
                        'end' => $end,
                        'arrow' => $arrow };
        } elsif ($shape eq "t") {
            $this = {'shape' => 't',
                     'ptx' => $x1,
                     'pty' => $y1,
                     'dir' => $arrow,
                     'color' => $palette->{$color}}
        } else { print "What is this? $shape\n" }

        push @{$MAP->[$depth]}, $this;
    }
    return $MAP;
}

sub Arrow {
    my $im = shift;
    my $ptx = shift;
    my $pty = shift;
    my $l = shift;
    my $w = shift;
    my $dir = shift;
    my $fill = shift;
    my $rad = rad($dir);
    my $rx = $ptx + int((sin($rad) * $l) - (cos($rad) * (0.5 * $l)));
    my $ry = $pty + int((cos($rad) * $l) + (sin($rad) * (0.5 * $l)));
    my $lx = $ptx + int((sin($rad) * $l) + (cos($rad) * (0.5 * $l)));
    my $ly = $pty + int((cos($rad) * $l) - (sin($rad) * (0.5 * $l)));
    #print "Arrow: $ptx/$pty; $rx/$ry; $lx/$ly  $dir $rad   $l/$w  $fill\n";

    my $tri = new GD::Polygon;
    $tri->addPt($ptx, $pty);
    $tri->addPt($rx, $ry);
    $tri->addPt($lx, $ly);

    $im->filledPolygon($tri, $fill);
}

1;