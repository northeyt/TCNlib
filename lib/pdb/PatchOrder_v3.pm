package PatchOrder_v3;

use strict;
use warnings;
use Math::Trig;

use Exporter;

our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

$VERSION = 1.00;
@ISA     = qw(Exporter);
@EXPORT  = ();
@EXPORT_OK = qw(PatchOrder);

use Carp;

sub PatchOrder {

    my($v, $patch, @atomLines) = @_;
    
    my ($centralRes, @patchResNumbers)
        = map { $_ =~ /[A-Z]:*(-*\d+[A-Z]*)/ig }(split(' ', $patch));
    # Capture negative res ids and A-Z insertion codes            

    my ($chain) = $patch =~ /patch ([A-Z])/; 
    
    print "The central residue is $centralRes\n" if $v;
            
    my %xyzHash = ();
    # Second hash is needed to later identify which amino acid corresponds
    # to which residue number
            
    my %aacidType = ();
    my $skip = 0; # Flag to skip patch if it contains unnatural amino acids
            
    # Get C-alpha co-ordinates
    foreach my $residue (@patchResNumbers) {
        foreach my $line (@atomLines) { 
            if ($line =~ /ATOM\s*\d+\s*CA\s+\w+\s*$chain\s*$residue\s+/) {
                                
                print "$residue: $line" if $v;
                
                my $x = substr($line, 31 - 1, 8);
                $x =~ s/\s+//g;
                my $y = substr($line, 39 - 1, 8);
                $y =~ s/\s+//g;
                my $z = substr($line, 47 - 1, 8);
                $z =~ s/\s+//g;
                $xyzHash{$residue} = [$x, $y, $z];
                                
                print "$x, $y, $z\n" if $v;
                
                # Get 1 Letter Code for the amino acid type (used later)
                $aacidType{$residue} = CodeTransform(substr($line, 17, 3));
                unless ($aacidType{$residue}) {
                    $skip = 1;
                }
            }
        }
    }
    
    if ($skip) {
        croak "Non-natural amino acid found in patch";
    }
    
    if (@patchResNumbers == 1) { # If patch only includes central residue
        return $aacidType{$patchResNumbers[0]};
    }
    elsif (@patchResNumbers == 2) {
        my ($neCentral) = grep(!/^$centralRes$/, @patchResNumbers);    
        return  "$aacidType{$centralRes}$aacidType{$neCentral}";
    }
    elsif (@patchResNumbers == 3) {
        my(@neCentral) = grep(!/^$centralRes$/, @patchResNumbers);
        for(my $i = 0; $i < @neCentral ; ++$i){
            $neCentral[$i] = $aacidType{$neCentral[$i]};
        }
        my $standardRep = StandardRep(join('', @neCentral));
        return "$aacidType{$centralRes}$standardRep";
    }
            
    # Central residue is i
    my($xi, $yi, $zi) = @{$xyzHash{$centralRes}};
            
    # subtract centralRes(x y z) from all patch residues
    #( move patch to origin), calculate magnitudes and
    #f ind the most distant point from origin.
    
            print "Here are the points translated so that $centralRes is at the origin\n" if $v;
            
    foreach(keys %xyzHash){
        
        print "BEFORE: $_ => @{$xyzHash{$_}}\n" if $v;
        
        # Subtract centralRes(x y z) from residue(x y z)
        @{ $xyzHash{$_} }
            = ( $xyzHash{$_}->[0] - $xi,
                $xyzHash{$_}->[1] - $yi,
                $xyzHash{$_}->[2] - $zi);
                    
        # To the end of each array, add the vector magnitude
        push( @{ $xyzHash{$_} },
              sqrt(  $xyzHash{$_}->[0]**2
                   + $xyzHash{$_}->[1]**2
                   + $xyzHash{$_}->[2]**2 )
          );
                    
        print "AFTER: $_ => @{$xyzHash{$_}}\n" if $v;
    }
            
    my @sorted = sort { $xyzHash{$b}->[3] <=> $xyzHash{$a}->[3] }
        (keys %xyzHash);
    
    my $primMag = $xyzHash{$sorted[0]}->[3];
    my $primRes = $sorted[0];
    my $scndMag =  $xyzHash{$sorted[1]}->[3];
    my $scndRes =  $sorted[1];
            
    print "Prim. reference residue is $primRes, magnitude $primMag\n"
        if $v;
    print "Second.  reference residue is $scndRes, magnitude $scndMag\n"
        if $v;
            
            
    # Get rotation angles for refRes
    
    #Get angle for rotation about x
    my($zyAngle) = TrueRotationAngle($xyzHash{$primRes}->[1],
                                     $xyzHash{$primRes}->[2]);
            
    print "The first rotation angle about x is $zyAngle\n" if $v;
    
    # Calculate rotation matrices for rotation about x
    my($rotX) = RotationMatrices($zyAngle, 'x');
            
    print "Here is the rotation about x matrix: @$rotX\n" if $v;
            
    # Use rotation matrices to transform y and z co-ordinates of all  points
    foreach (keys %xyzHash) {
        
        print "Original co-ordinates for $_ : @{$xyzHash{$_}} \n" if $v;
        
        @{ $xyzHash{$_} }
            = DoRotation( $xyzHash{$_}->[0],
                          $xyzHash{$_}->[1],
                          $xyzHash{$_}->[2],
                          $rotX);
                }
            
    # Now calculate angle for rotation about y
    my($xzAngle) = TrueRotationAngle($xyzHash{$primRes}->[2],
                                     $xyzHash{$primRes}->[0]);
            
    # Calculate rotation matrices for rotation about y
    my($rotY) = RotationMatrices($xzAngle, 'y');
    
    # Use rotation matrices to transform x and z co-ordinates of all points
    foreach (keys %xyzHash) {
        print "New co-ordinates for $_ : @{$xyzHash{$_}} \n" if $v;
        
        @{ $xyzHash{$_} }
            = DoRotation($xyzHash{$_}->[0],
                         $xyzHash{$_}->[1],
                         $xyzHash{$_}->[2], $rotY);
        
        print "Final co-ordinates for $_ : @{$xyzHash{$_}} \n" if $v;
    }
            
    # Perform a third rotation, about the x axis,
    # using the second reference residue
    
    $zyAngle = TrueRotationAngle($xyzHash{$scndRes}->[1],
                                 $xyzHash{$scndRes}->[2]);
            
    # Calculate rotation matrices for rotation about x
    $rotX = RotationMatrices($zyAngle, 'x');
            
    # Use rotation matrices to transform x and z co-ordinates of all points
    foreach (keys %xyzHash) {
        print "New co-ordinates for $_ : @{$xyzHash{$_}} \n" if $v;
        @{$xyzHash{$_}}
            = DoRotation($xyzHash{$_}->[0],
                         $xyzHash{$_}->[1],
                         $xyzHash{$_}->[2],
                         $rotX);
        
        print "Final co-ordinates for $_ : @{$xyzHash{$_}} \n" if $v;
    }
            
    #Initialize unit vectors
            
    # Unit vector h, in direction of x axis
    my($hx, $hy, $hz) = (1, 0, 0);
    
    # Magnitude of unit vector h
    my $absh = 1;
            
    # Unit vector i, in direction of y axis
    my ($ix, $iy, $iz) = (0, 1, 0);
    
    # Magnitude of unit vector i
    my $absi = 1;
    
    # Remove central residue from hash before calculating angles
    delete $xyzHash{$centralRes};
    
    my %angleHash = ();
    
    foreach  (keys %xyzHash) {
        print "Residue $_ \n" if $v;
        
        my($kx, $ky, $kz) = @{$xyzHash{$_}};
        
        # Calculate k vector magnitude
        my $absk = sqrt($kx * $kx + $ky * $ky);
        
        # Calculate cos Angle using inner-product and dot products of h and k
        my $cosAngle = ($hx * $kx + $hy * $ky) / ($absh * $absk);
        print "The cosAngle is $cosAngle\n" if $v;
                    
        # Calculate cosAngleCheck using inner-product and
        # dot products i and k
        
        my $cosAngleCheck = ($ix * $kx + $iy * $ky) / ($absi * $absk);
        print "The cosAngleCheck is $cosAngleCheck\n" if $v;
                    
        # Compare the signs of cosAngle and cosAngleCheck to determine
        # reflex angles
                    
        my $angle = 0;
        
        if($cosAngle > 0) {
            if($cosAngleCheck <  0) {
                $angle = (2 * pi)  - acos($cosAngle);
            }
            else {    
                $angle = acos($cosAngle);
            }
        }
        else{
            
            print "$_ cosAngle = $cosAngle, cosAngleCheck = $cosAngleCheck\n" if $v;
            
            if($cosAngleCheck < 0) {
                $angle = (2 * pi) - acos($cosAngle);
            }
            else {
                $angle = acos($cosAngle);
            }
        }
        
        $angleHash{$_} = $angle;
    }
            
    # Now sort the residues by angle to determine the order around
    # the central residue
    
    # Then, make an array with the central residue as the first element,
    # and the ordered amino acids as the second element
            
    my @patchArray = ();
    $patchArray[0] = $aacidType{$centralRes};

    my @sorted_residues
        = sort { $angleHash{$a} <=> $angleHash{$b} } (keys %angleHash); 
    
    foreach my $residue (@sorted_residues){
        print "$residue($aacidType{$residue}) " if $v;
        $patchArray[1] .= $aacidType{$residue};
    }
    
    print "\n@patchArray \n" if $v;
    
    # Now standardize the second element for comparison to other
    # patch descriptors.
    
    $patchArray[1] = StandardRep($patchArray[1]);
    
    print "@patchArray \n" if $v;
    
    my $order = "$patchArray[0]$patchArray[1]";
    $order =~ s/\s//g;

    return $order;
}
    


sub DoRotation
{

    my($x, $y, $z, $rMatrixRef) = @_;

    my @rMat = @{$rMatrixRef};

    my $xd = ($x*$rMat[0] + $y*$rMat[1] + $z*$rMat[2]);
    my $yd = ($x*$rMat[3] + $y*$rMat[4] + $z*$rMat[5]);
    my $zd = ($x*$rMat[6] + $y*$rMat[7] + $z*$rMat[8]);

    return($xd, $yd, $zd);

}

sub RotationMatrices # Give this subroutine the angle and the axis that you want to rotate about
{

    my($angle, $axis) = @_;
    $axis = lc($axis);
    my @rotationXMatrix = (0);
    my @rotationYMatrix = (0);
    my @rotationZMatrix = (0);

    if($axis eq 'x')
    {
	@rotationXMatrix =
	(
	 1, 0, 0,
	 0, cos($angle), -sin($angle),
	 0, sin($angle), cos($angle)
	);

	return \@rotationXMatrix;

    }
    elsif($axis eq 'y')
    {
	@rotationYMatrix =
	(
	 cos($angle), 0, sin($angle),
	 0, 1, 0,
	 -sin($angle), 0, cos($angle)
	);

	return \@rotationYMatrix;

    }
    elsif($axis eq 'z')
    {
	@rotationZMatrix =
	(
	 cos($angle), -sin($angle), 0,
	 sin($angle), cos($angle), 0,
	 0, 0, 1
	);

	return \@rotationZMatrix;

    }
    else
    {
	die print "ERROR: Axis for rotation is not valid!\n";
    }
}


# Calculates atan(opp/adj) and real angle value according to opp and adj sign.
sub TrueRotationAngle
{

    my($opp, $adj) = @_;
    use Math::Trig;

    my $angle = 0;

    if ($adj != 0.0)
    {
	$angle = atan($opp / $adj);

	if($opp < 0.0 && $adj > 0.0)
	{
	    $angle += 2*pi;
	}
	if($adj < 0.0)
	{
	    $angle += pi;
	}
    }
    else
    {
	if($opp > 0.0) # #1st->2nd quadrant boundary
	{
	    $angle = pi / 2;
	}
	else  # 3rd -> 4th quadrant boundary
	{
	    $angle = 3 * (pi / 2);
	}
    }

    if($opp == 0)
    {
	if($adj > 0.0) # 4th -> 1st quadrant boundary
	{
	    $angle = 0;
	}
	else # 2nd -> 3rd quadrant boundary
	{
	    $angle = pi;
	}
    }

    return($angle);

}

sub CodeTransform
{
    my ($acid) = @_;

    $acid = uc($acid);

    my %acidHash = (
		    'ALA' => 'A',
		    'VAL' => 'V',
		    'LEU' => 'L',
		    'ILE' => 'I',
		    'PRO' => 'P',
		    'TRP' => 'W',
		    'PHE' => 'F',
		    'MET' => 'M',
		    'GLY' => 'G',
		    'SER' => 'S',
		    'THR' => 'T',
		    'TYR' => 'Y',
		    'CYS' => 'C',
		    'ASN' => 'N',
		    'GLN' => 'Q',
		    'LYS' => 'K',
		    'ARG' => 'R',
		    'HIS' => 'H',
		    'ASP' => 'D',
		    'GLU' => 'E',
		   );

    if (defined $acidHash{$acid}) {
        return $acidHash{$acid};
    }else{
        return 0;
    }
}

sub StandardRep
    {
        my ($residues) = @_;
        
        my $standardRep = "Z";

        # Double up residue string for comparison
        my $dResidues = $residues.$residues;

        for (my $i = 0 ; $i < (length $residues) ; $i++)
            {
                if(substr($dResidues, $i, length $residues) lt $standardRep)
                    {
                        $standardRep = substr($dResidues, $i, length $residues);
                    }
            }
        return $standardRep;
    }

sub Test
    {

        my %xyzHash = @_;
        
        
        # TEST: perform series of x y and z  rotations on all points to see if this has an effect on order
        
        # Calculate rotation matrices for rotation about x
        
        my($rotationTestMatrixref) = RotationMatrices((0), 'x');
        
        # Use rotation matrices to transform co-ordinates of all points
        
        foreach (keys %xyzHash)
            
            {
                print "Here are the new co-ordinates for $_ : @{$xyzHash{$_}} \n";
                @{$xyzHash{$_}} = DoRotation(@{$xyzHash{$_}}[0], @{$xyzHash{$_}}[1], @{$xyzHash{$_}}[2], $rotationTestMatrixref);
                print "Here are the final co-ordinates for $_ : @{$xyzHash{$_}} \n";
            }
        
        # Calculate rotation matrices for rotation about y
        
        $rotationTestMatrixref = RotationMatrices((0), 'y');
        
        # Use rotation matrices to transform co-ordinates of all points
        
        foreach (keys %xyzHash)
            
            {
                print "Here are the new co-ordinates for $_ : @{$xyzHash{$_}} \n";
                @{$xyzHash{$_}} = DoRotation(@{$xyzHash{$_}}[0], @{$xyzHash{$_}}[1], @{$xyzHash{$_}}[2], $rotationTestMatrixref);
                print "Here are the final co-ordinates for $_ : @{$xyzHash{$_}} \n";
            }
        
        # Calculate rotation matrices for rotation about z
        
        $rotationTestMatrixref = RotationMatrices((0), 'z');
        
        # Use rotation matrices to transform co-ordinates of all points
        
        foreach (keys %xyzHash)
            
            {
                print "Here are the new co-ordinates for $_ : @{$xyzHash{$_}} \n";
                @{$xyzHash{$_}} = DoRotation(@{$xyzHash{$_}}[0], @{$xyzHash{$_}}[1], @{$xyzHash{$_}}[2], $rotationTestMatrixref);
                print "Here are the final co-ordinates for $_ : @{$xyzHash{$_}} \n";
            }
    }
   
sub Residue2Group {
    my($residue) = @_;

    # A = hydrophobic, B = polar non-charged, C = positive, D = negative
    
    my %hash =
        (
            'G' => 'H',
            'A' => 'H',
            'V' => 'H',
            'L' => 'H',
            'I' => 'H',
            'M' => 'H',
            'F' => 'H',
            'W' => 'H',
            'P' => 'H',
            'S' => 'B',
            'T' => 'B',
            'C' => 'B',
            'Q' => 'B',
            'N' => 'B',
            'Y' => 'B',
            'D' => 'D',
            'E' => 'D',
            'K' => 'C',
            'R' => 'C',
            'H' => 'C'
        );
    return $hash{$residue};
}

1;
