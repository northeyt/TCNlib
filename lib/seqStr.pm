package seqStr;

sub parse {
    my $seqStr = shift;
    if ($seqStr =~ /^([A-Z-\n\s]+)$/im){
        my $cleanString = $1;
        $cleanString =~ s/\s//gms;
        return $cleanString;
    }
    else {
        return;
    }
}

sub isFASTA {
    return $_[0] =~ />(.*)\n([A-Z-\n\s]+)/i;
};

sub parseSeqFromFASTA {
    my $str = shift;
    if (my ($id, $seqString) = isFASTA($str)){
        return parse($seqString);
    }
    else {
        return;
    }
}

sub parseIDFromFASTA {
    my $str = shift;
    if (my ($id) = isFASTA($str)){
        return $id;
    }
    else {
        return;
    }
}

1;
