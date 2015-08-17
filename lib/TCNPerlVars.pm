package TCNPerlVars;
# Where this file lives
$libdir = "~/scripts/lib/tcnutils";

# Data
$datadir = "/acrm/data/people/zcbtfo4";

## bioplib and bioptools
$bioptoolsBin = $ENV{HOME} ."/bin";
$bioplibDataDir = $ENV{HOME} . "/data";

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

# pdb file cache
$pdb_file_cache_dir = $ENV{HOME} . '/pdb-cache';

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

# Command-line program paths for pdb and related classes
$makepatch     = "/home/bsm/martin/bin/makepatch";
$xmas2pdb      = "/home/bsm/martin/acrm/CONSULTANCY/inpharmatica/software/xmas2pdb/xmas2pdb";
$pdb2xmas      = "/home/bsm/martin/acrm/CONSULTANCY/inpharmatica/software/pdb2xmas/pdb2xmas";
$pdb2xmas_bin  = "/home/bsm/martin/acrm/CONSULTANCY/inpharmatica/software/bin/";
$idabchain     = "$bioptoolsBin/idabchain";
$kabatnum      = '/home/bsm/martin/abnum/installed/numbering/kabatnum.pl';
$chaincontacts = "$bioptoolsBin/chaincontacts";
$blastall      = "/Users/tcn/software/blast-2.2.22/bin/blastall";
$princip64     = "/acrm/usr/local/apps/surfnet/bin64/princip";
$scorecons     = "$bioptoolsBin/scorecons";
$pdbsolv       = "$bioptoolsBin/pdbsolv";
$pdbsslist     = "$bioptoolsBin/pdblistss";
$pdbsecstr     = "$bioptoolsBin/pdbsecstr";
$pdbhbond      = "$bioptoolsBin/pdbhbond";

# Data files for pdb and related classes
$hydroPhoValueFile = $bioplibDataDir . "/kyte.hpb";
$radii_file = "$bioplibDataDir/radii.dat";
$asurf_radii_file = "/acrm/data/people/zcbtfo4/vdw.radii";

# CDhit and related
$cdhit = '/home/bsm3/zcbtfo4/cd-hit-4.6.1/cd-hit';
$clstr2xml = '/home/bsm3/zcbtfo4/cd-hit-4.6.1/clstr2xml.pl';

# Clustal and Muscle (MSAs)
$clustalw = '/acrm/usr/local/bin/clustalw';
$clustalO = '/home/bsm3/zcbtfo4/clustalo-1.2.0-Ubuntu-x86_64';
$muscle   = $ENV{'HOME'} . '/software/muscle-3.8.31/muscle';

# asurf64 and related
$asurf64 = "/home/bsm/martin/bin/asurf64";
$standardData = "$libdir/pdb/standard.data";

# SACS Antibody-containing PDB XML File 
$SACSxml = '/acrm/data/abs/xml/all.xml';
    
# My normal binaries directory
$bindir  = "/home/bsm/martin/bin";

# Other binaries directories
$ssapbindir = "/acrm/home/andrew/ssap/bin";
$mlsabindir = "/acrm/home/andrew/mlsa";
$sstrucbindir = "/acrm/home/andrew/sstruc/bin";

# Environment variables for Kabat related programs 
$ENV{'KABATALIGN'} = $bioplibDataDir; # Alignment matrices

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

# WEKA related
$wekaLib     = $ENV{'HOME'} . "/software/weka-3-7-10/weka.jar";
$javaForWeka = "/usr/bin/java";

# seq dbs formatted for blast searching
$swissProtDB  = $ENV{'HOME'} . '/data/blast/swissprot';
$pdb_db       = $ENV{'HOME'} . '/data/blast/pdbaa';

1;

