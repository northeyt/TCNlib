package TCNPerlVars;

# Where this file lives
$libdir = "/home/bsm3/zcbtfo4/scripts/lib";

# Data
$datadir = "/acrm/data/people/zcbtfo4";

# Tmp dir
$tmpdir = "/tmp/TCN";

# About the PDB
$pdbdir       = "/acrm/data/pdb";
$pdbprep      = "/acrm/data/pdb/pdb";
$pdbprepname  = "pdb";
$pdbext       = ".ent";

# Obsolete PDB
$obspdbdir = "/acrm/data/pdb_obsolete/uncompressed_files";
$obspdbprep = "pdb";
$obspdbext = ".ent";

# Other PDB-related files
$pdbdomdir = "/acrm/data/dompdb";
$domsstdir = "/acrm/data/domsst";
$sstdir    = "/acrm/data/sst";
$pdbseqdir = "/acrm/data/pdbseq";
$pisadir   = "/acrm/data/pisa";
$pisaext   = ".pisa";
$xmasdir   = "/acrm/data/xmas/pdb";
$xmasprep  = "/acrm/data/xmas/pdb/pdb";
$xmasext   = ".xmas";

# makepatch and related
$makepatch = "/home/bsm/martin/bin/makepatch";
$xmas2pdb  = "/home/bsm/martin/acrm/CONSULTANCY/inpharmatica/software/xmas2pdb/xmas2pdb";

$radii_file = "/acrm/data/people/zcbtfo4/radii.dat";
$asurf_radii_file = "/acrm/data/people/zcbtfo4/vdw.radii";
$pdb2xmas   = "/home/bsm/martin/acrm/CONSULTANCY/inpharmatica/software/pdb2xmas/pdb2xmas";

$pdb2xmas_bin
    = "/home/bsm/martin/acrm/CONSULTANCY/inpharmatica/software/bin/";

# Command-line program paths for pdb and related classes
$getresol = '/home/bsm/martin/bin/getresol';
$idabchain = '/home/bsm/martin/bin/idabchain';
$kabatnum = '/home/bsm/martin/abnum/installed/numbering/kabatnum.pl';
$chaincontacts = '/home/bsm/martin/bin/chaincontacts';
$pdb_db = "/home/bsm3/zcbtfo4/data/pdb/db/pdbaa";
$blastall = "/acrm/usr/local/bin/blastall";

# CDhit and related
$cdhit = '/home/bsm3/zcbtfo4/cd-hit-4.6.1/cd-hit';
$clstr2xml = '/home/bsm3/zcbtfo4/cd-hit-4.6.1/clstr2xml.pl';

# ClustalW
$clustalw = '/acrm/usr/local/bin/clustalw';
$clustalO = '/home/bsm3/zcbtfo4/clustalo-1.2.0-Ubuntu-x86_64';

# asurf64 and related
$asurf64 = "/home/bsm/martin/bin/asurf64";
$standardData = "/home/bsm3/zcbtfo4/scripts/lib/pdb/standard.data";

# SACS Antibody-containing PDB XML File 
$SACSxml = '/acrm/data/abs/xml/all.xml';
    
# My normal binaries directory
$bindir  = "/home/bsm/martin/bin";

# Other binaries directories
$ssapbindir = "/acrm/home/andrew/ssap/bin";
$mlsabindir = "/acrm/home/andrew/mlsa";
$sstrucbindir = "/acrm/home/andrew/sstruc/bin";

# Environment variables for Kabat related programs
#$ENV{'KABATALIGN'} = "/acrm/home/andrew/kabat/data"; # Alignment matrices
$ENV{'KABATDIR'}   = "/acrm/data/kabat/kabatman";    # KabatMan data 
$ENV{'KABATALIGN'} = "/home/bsm/martin/kabat/data"; # Alignment matrices
#$ENV{'KABATDIR'}   = "/home/bsm/martin/kabat/data";    # KabatMan data 



# My general data directory
$ENV{'DATADIR'}    = "/home/bsm/martin/data";


# BLAST related
$BlastAllExe       = "/acrm/usr/local/bin/blastall";
$BlastPGPExe       = "/acrm/usr/local/bin/blastpgp";
$BlastDBDir        = "/acrm/data/blastdb";
$BlastRCFile       = "/acrm/home/andrew/.ncbirc";

# FASTA related
$ssearch = "/acrm/usr/local/bin/ssearch33";
$fasta   = "/acrm/usr/local/bin/fasta33";
$ssearch_64 = "/acrm/usr64/local/bin/ssearch33";
$fasta_64   = "/acrm/usr64/local/bin/fasta33";

# PostgreSQL related
$pghost = "acrm8";
$psql   = "/acrm/usr/local/bin/psql";

$ENV{'PGHOST'} = $pghost;
$ENV{'PGLIB'} = "/usr/lib/pgsql";
$ENV{'LD_LIBRARY_PATH'} = "$ENV{'LD_LIBRARY_PATH'}:/usr/lib/pgsql";

# SAAP related
$saapServerBindir = "/home/bsm/martin/SAAP/server/";


1;

