package pdb::ViewPatch;
use Moose;
use Carp;
use TCNUtil::types;
use TCNUtil::ARFF;
use TCNUtil::write2tmp;
use pdb;

has 'pymolPipe' => (
    is => 'rw',
    isa => 'GlobRef',
    lazy => 1,
    builder => '_getPymolPipe'
);

has 'chainObjectLookup' => (
    is => 'rw',
    isa => 'HashRef',
    default => sub {{}},
);

has 'patchDir' => (
    is => 'rw',
    isa => 'Directory',
    predicate => 'has_patchDir',
);

has 'patchDescriptionLookup' => (
    is => 'rw',
    isa => 'HashRef',
    lazy => 1,
    builder => '_buildLookupFromPatchDir',
);

has 'arff' => (
    is => 'rw',
    isa => 'ARFF',
    predicate => 'has_arff',
);

sub view {
    my $self = shift;
    my $patchID = shift;
    my ($pdbCode, $chainID, $centralResID) = $self->_parsePatchID($patchID);
    my $chainObjectID = $self->prepareChainObject("$pdbCode:$chainID");
    $self->_prepareChainObjectDisplay($chainObjectID);
    my $patchSelectionID = join("_", ($pdbCode, $chainID, $centralResID));
    my $createPatchObjectCmd
        = $self->_prepareSelectPatchCommand($chainObjectID,
                                            $patchID,
                                            $patchSelectionID);
    $self->runCmd($createPatchObjectCmd);
    $self->_preparePatchObjectDisplay($patchSelectionID);
    if ($self->has_arff) {
        my $instance = $self->arff->findInstance(sub {$_->id eq $patchID});
        if (! $instance) {
            print "No instance found in arff with id $patchID\n";
        }
        else {
            my %attrs = $instance->getHashOfAttributeNameToValue();
            print "Patch $patchID\n";
            $self->printAttributes(%attrs);
        }
    }
}

sub printAttributes {
    my $self = shift;
    my %attrs = @_;
    print "$_ = " . $attrs{$_} . "\n" foreach grep {$_ ne "patchID"} sort {$a cmp $b} keys %attrs;
}


sub _buildLookupFromPatchDir {
    my $self = shift;
    croak "No patchDir has been set!" if ! $self->has_patchDir;
    opendir(my $DIR, $self->patchDir)
        or die "Cannot open patch dir " . $self->patchDir . " $!";
    my %patchID2Summary = ();
    
    while (my $fName = readdir($DIR)) {
        next if $fName =~ /^\./; # Skip . and ..
        my $file = $self->patchDir . "/$fName";
        open(my $IN, "<", $file) or die "Cannot open file $file, $!";
        my ($pdbID) = $fName =~ /(\S+)\.patches/g;
        my $pdbCode = substr($pdbID, 0, 4);
        foreach my $summaryLine (<$IN>) {
            chomp $summaryLine;
            my ($cResID, @resIDs) = patch::parseSummaryLine($summaryLine);
            unshift(@resIDs, $cResID);
            $cResID =~ s/\./:/;
            $patchID2Summary{"$pdbCode:$cResID"} = \@resIDs;
        } 
    }
    return \%patchID2Summary;
}

sub _parsePatchID {
    my $self = shift;
    my $patchID = shift;
    my ($pdbCode, $chainID, $cResSeq) = split(/:/, $patchID);
    return ($pdbCode, $chainID, $cResSeq);
}

sub _getPymolPipe {
    my $pymol = `which pymol`;
    chomp $pymol;
    open(my $OUT, "|$pymol -p > /dev/null") or die "Pipe failed\n";
    return $OUT;
}

sub _prepareChainObjectDisplay {
    my $self = shift;
    my $chainObjectID = shift;
    $self->runCmd("hide everything");
    $self->runCmd("show surface, $chainObjectID");
    $self->runCmd("color red, $chainObjectID");
}

sub _preparePatchObjectDisplay {
    my $self = shift;
    my $patchSelectionID = shift;
    $self->runCmd("show surface, $patchSelectionID");
    $self->runCmd("color blue, $patchSelectionID");
}

sub _prepareSelectPatchCommand {
    my $self = shift;
    my $chainObjectID = shift;
    my $patchID = shift;
    my $patchSelectionID = shift;
    my $resSelection  = $self->_getSelectionForPatchID($patchID);
    return "select $patchSelectionID, ($chainObjectID and ($resSelection))";
}

sub runCmd {
    my $self = shift;
    my $cmd  = shift;
    {
        my $oldOut = select;
        select($self->pymolPipe);
        local $| = 1;
        print "$cmd\n";
        select($oldOut);
    }
}

sub _getSelectionForPatchID {
    my $self = shift;
    my $patchID = shift;
    croak {message => "There is no patch summary for patch $patchID!"}
        if ! exists $self->patchDescriptionLookup->{$patchID};
    
    my @resSelections = map {$self->_selectionStringFromResID($_)}
        @{$self->patchDescriptionLookup->{$patchID}};
    my $selectionString = join (" or ", @resSelections);
    return $selectionString;
}

sub _selectionStringFromResID {
    my $self  = shift;
    my $resID = shift;
    my ($chainID, $resSeq) = split(/\./, $resID);
    return "(chain $chainID and resi $resSeq)";
}

sub prepareChainObject {
    my $self = shift;
    my $pdbID = shift;
    if (! exists $self->chainObjectLookup->{$pdbID}) {
        my ($pdbCode, $chainID) = split(/:/, $pdbID);
        my $chain = chain->new(pdb_code => $pdbCode, chain_id => $chainID);
        my $file = write2tmp->new(data => [map {"$_"} @{$chain->atom_array}])->file_name();
        my $objectID = $pdbID;
        $objectID =~ s/:/_/;
        $self->runCmd("load $file, $objectID");
        $self->chainObjectLookup->{$pdbID} = $objectID;
    }
    return $self->chainObjectLookup->{$pdbID};
}

1;
