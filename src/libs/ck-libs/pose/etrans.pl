#!/usr/bin/perl
#!/bin/sh 
#exec perl -w -x $0 ${1+"$@"} # -*- mode: perl; perl-indent-level: 2; -*-
#!perl -w 

#$Id$

use FileHandle;
use English;
use strict; 

my $infile = $ARGV[0];
my $inci = "$infile.ci";
my $inh = "$infile.h";
my $inC = "$infile.C";
my $outci = "$infile\_sim.ci";
my $outh = "$infile\_sim.h";
my $outC = "$infile\_sim.C";
my (@strategies,@representations,@posera,@groups,@newline);
my (%posers,%strat, %methods,%rep,%base,%events,%nonevents,%messages,%needsopeq,%puppy,%cppuppy,%emessages);
my ($class,$method,$message);
my ($template,$tdef, $classpattern); #for template support
my $wantindent=0;
print "Processing .ci files...\n";
#hashes for parse state
my %openblock;
my %closeblock;
my %closecomment;
my %opencomment;
my %incomment;
my %blocklevel;
my %lastline;
# BEGIN READ HEADERS & MESSAGES
#read inci
my $incihandle= new FileHandle();
$incihandle->open("$inci") or die "cannot open $inci";
inithandle($incihandle);
my $outcihandle=new FileHandle();
$outcihandle->open(">$outci") or die "cannot open $outci";
my @line;
my $thisline;
my @otherstuff;
while (@line=split(' ',($thisline=getcodeline($incihandle)))) {
  if ($line[0] eq "message") {
    chop $line[1];
    $messages{$line[1]}=$lastline{$incihandle};
  }
  elsif ($line[0] eq "poser") {
    $posers{$line[1]}=poserdeal($incihandle,\@line,$thisline);
    push(@posera,$posers{$line[1]});
  }
  elsif ($line[0] =~ /((group)|(chare))/) {
    push(@groups,groupnchardeal($incihandle,\@line,$thisline));
  }
  else
  {
      push(@otherstuff,$thisline);
  }
}

foreach my $j (keys %events) {
  foreach my $i (@{$events{$j}}) {
    print "class $j has method $i->[0] $i->[1] \n";
  }
}


$incihandle->close();
%emessages=%messages;
# on to classes
# END READ HEADERS & MESSAGES

# BEGIN HEADERS
# write outci
print "Processing headers...\n";
$outcihandle->print("module $infile {\n");
$outcihandle->print("  extern module sim;\n");
$outcihandle->print("  extern module pose;\n");
$outcihandle->print(@otherstuff);
$outcihandle->print("\n");
foreach my $i (keys %messages) {
  $outcihandle->print("  message $i;\n");
}
$outcihandle->print("\n");

# write outh
my $outhhandle=new FileHandle();
my $indent=`indent --version`;
if ( $wantindent && ( $indent =~ /GNU/) && ( $indent !~ /2.2.5/) ) {
  $outhhandle->open("|indent -npcs -npsl -bap -bad -br -nce > $outh") or die "cannot open indented $outh";
} else {
  $outhhandle->open("> $outh") or die "cannot open $outh";
}
$outhhandle->print("#ifndef $infile\_H\n");
$outhhandle->print( "#define $infile\_H\n");
$outhhandle->print( "#include \"pose.h\"\n");
$outhhandle->print( "#include \"$infile.decl.h\"\n");
$outhhandle->print( "\n");
#foreach my $i (@strategies) {
#    $outhhandle->print( "#include \"$i.h\"\n");
#}
#foreach my $i (@representations) {
#    if ($i ne "rep") {
#	$outhhandle->print( "#include \"$i.h\"\n");
#    }
#}

# END HEADERS

# BEGIN READ CLASSES
print "Processing classes...\n";
my $group;
foreach $group (@groups) {
  foreach (@{$group->{literal}}) {
    $outcihandle->print($_."\n");
  }
}
my $poser;
foreach $poser (@posera) {
  foreach (@{$poser->{literal}}) {
    $outcihandle->print($_."\n");
  }
}
#footer
$outcihandle->print("};\n"); 
$outcihandle->close;

# PRESCAN $infile.h for cpPups.  Record them by class so we don't make new ones.
my $inhhandle=new FileHandle();
$inhhandle->open("$inh") or die "cannot open $inh";
inithandle($inhhandle);
while (@line=split(' ',($thisline=getcodeline($inhhandle))))
  {
    #seek the puppy
    if ($thisline =~ /^\s*class\s+([^:\s\(\{]+)/)
      {
	$class=$1;
      }
    if ($thisline =~ /void\s+cpPup\s*\(PUP\s*/)
      {
	$cppuppy{$class}=1;
	print "found user defined cppup for class $class \n";
      }
  }
# BEGIN READ MESSAGES FROM $infile.H
print "Processing messages...\n";
$inhhandle->open("$inh") or die "cannot open $inh";
inithandle($inhhandle);
my ($currep,$curstrat,$annoying,$declaration);
my ($opeq,$destructor);
my ($numClasses);
my (@pair);
# foreach message append :public eventMsg to decl
# reject message if no operator= defined within
# foreach class save body in array and make two copies
# first one is written verbatim
# second one append :public sim to header decl
my $inbody=0;
my $inpup=0;
my ($ismessage,$issim,$isconstructor);
my @body;
while (@line=split(' ',($thisline=getcodeline($inhhandle)))) {
  if ($inbody) {
    if (posefunctest($thisline)>0) {
      push(@body,posefuncmap($thisline,$issim,$isconstructor,$class));
    } else {
      push(@body,$lastline{$inhhandle});
    }
    if ($thisline =~ /$class&*\s*operator\s*=.*$class&/i) {
      $opeq=1;
    }
    if ($thisline =~ /~$class/i) {
      $destructor=1;
    }
    #seek the puppy
    if ($thisline =~ /void\s+pup\s*\(PUP\s*::\s*er/) {
      $puppy{$class}=[];
      #	  print "Found puppy for $class\n";
      #got our puppy
      if ($openblock{$inhhandle}) {
	# got the body here
	$inpup=$blocklevel{$inhhandle};
      } elsif ($currep eq "chpt") {
	push(@body,"  void cpPup(PUP::er &p);") if(!exists ($cppuppy{$class}) );
      } else {			#do nothing
      }
    } elsif (($inpup)&&($thisline !~/((chpt<state)|(chpt::pup))/)) {
      #	    print "Inserting pup ".$lastline{$inhhandle}." for $class \n";
      push(@{$puppy{$class}},$lastline{$inhhandle});
    }
  }
  if ($inpup && ($blocklevel{$inhhandle}<$inpup) && ($closeblock{$inhhandle}==1)) { #got the whole pup
    #add the cpPup
    if (!exists($cppuppy{$class}) && ($currep eq "chpt")) {
      $inpup=0;
      push(@body,"  void cpPup(PUP::er &p) {");
      push(@body,@{$puppy{$class}});
    }
  }
  if (($blocklevel{$inhhandle}==0) && ($closeblock{$inhhandle}==1)) { #got the whole thing
    print "class $class issim $issim ismessage $ismessage $inbody defined \n";
    $inbody=0;
    if ($ismessage) {
      $outhhandle->print("class $class : public eventMsg {\n");
      foreach (@body) {
	$outhhandle->print("$_\n");
      }
      $outhhandle->print("\n");
    } elsif ($issim) {
      #take care of first class copy
      $outhhandle->print("class $class : public ".$base{$class}." {\n");
      $outhhandle->print(" private:\n   void ResolveFn(int fnIdx, void *msg);\n") if ($#{$events{$class}}>=0);
      $outhhandle->print("   void ResolveCommitFn(int fnIdx, void *msg);\n") if($#{$events{$class}}>=0);
      $outhhandle->print(" public:\n");
      $outhhandle->print("   $class(CkMigrateMessage *) {};\n");
      $outhhandle->print("   void pup(PUP::er &p);\n") if exists $puppy{$class};
      #	  $outhhandle->print("   void cpPup(PUP::er &p);\n");
      # Get verbatim from .ci file for SOUT
      # we stored that in the $posers hash
      my $thischar=$posers{$class};
      #for classes with strategies go get them
      foreach (@{$thischar->{literal}}) {
	my ($lmethod,$lmessage);
	my @newline = split;
	next if (/array \[1D\]/); #skip classdef
	next if (/\}\;/);	#skip close
	if ($newline[0] ne "entry") {
	  $outhhandle->print("$_\n");
	} else {
	  s/entry\s//;
	  s/\[.+\]\s//;
	  $outhhandle->print("  $_\n");
	  for my $i (@newline) {
	    if ($i =~ /([\w_]+)\((\w+)/) {
	      $lmethod = $1;
	      $lmessage = $2;
	      last;
	    }
	  }
	  if(defined $lmethod and defined $lmessage)
	    {
	      my $pair = [$lmethod, $lmessage];
	      push(@{$methods{$class}}, $pair);
	      print "\n  method  ".$class."::".$lmethod."(".$lmessage.")\n";
	    }
	}
      }
      $outhhandle->print("};\n\n");
      if (($curstrat eq "opt") || ($currep eq "chpt")) {
	print "warning ! no ::operator= in class $class" if (exists($needsopeq{$class}) &&( $opeq!=1));
	print "warning ! no ~ destructor in $class" if $destructor!=1;
      }
      if (!exists($puppy{$class})) {
	print "warning ! $class has no puppy, you loser \n";
      }
      #2nd copy
      if ($currep eq "chpt") {
	$outhhandle->print("class state_$class : public $currep<state_$class> {\n");
      } else {
	$outhhandle->print("class state_$class : public $currep {\n");
      }
      $outhhandle->print("  friend class $class;\n");
      my $declout=0;
      foreach (@body) {
	#substitute state_$class for occurence of $class
	s/(\b$class\b)/state_$1/g;
	$outhhandle->print("$_");
	$outhhandle->print("\n");
	#look for the public declaration

	if ((/public:/)&&(!$declout)) {
	  $outhhandle->print("  state_$class(sim *p, strat *s) { parent = p; myStrat = s; }\n");
	  $declout=1;
	}
      }
      $outhhandle->print("\n");
      #make state_$class constructor

    } else {			#dunno what this is, just pass it
      foreach (@body) {
	$outhhandle->print("$_");
	$outhhandle->print("\n");
      }
    }
    $inpup=0;
    @body=();
    $opeq=0;
    $destructor=0;
  } elsif (($blocklevel{$inhhandle}==1) && ($openblock{$inhhandle}==1) &&(!$inbody)) { # decl line
    # is  it a message or a sim class?
    # is a chpt?
#    print "declaration line $thisline \n";
    $inbody=1;
    if ($thisline =~ /^\s*class\s+([^:\s\(\{]+)/) {
      $class=$1;
      $tdef='';
      $classpattern=$class;
      if ($thisline =~ /:/) {
	$outhhandle->print($thisline."\n");
	$issim=0;
	$isconstructor=0;
	$ismessage=0;
	delete($emessages{$class}) if exists($emessages{$class});
      } elsif (defined($messages{$class})) {
	$ismessage=1;
	$issim=0;
      } elsif (exists($base{$class})) {
	$isconstructor=1;
	$issim=1;
	$ismessage=0;
	$curstrat=$strat{$class};
	$currep = $rep{$class};
	print "$curstrat... $currep...$class \n";
      } else {			#some sort of local use class
	$isconstructor=1;
	$issim=$ismessage=0;
	push(@body,$thisline);
      }
    } elsif ($thisline =~ /^\s*template\s*<\s*class\s+(.*)\s*>\s*class\s+([^:\s\(\{]+)/) {
      $class=$2;
      $tdef=$1;
      $classpattern=$class.'\s*<.*>';
      if ($thisline =~ /:/) {
	$outhhandle->print($thisline."\n");
	$issim=0;
	$isconstructor=0;
	$ismessage=0;
	delete($emessages{$class}) if exists($emessages{$class});
      } elsif (defined($messages{$class})) {
	$ismessage=1;
	$issim=0;
      } elsif (exists($base{$class})) {
	$isconstructor=1;
	$issim=1;
	$ismessage=0;
	$curstrat=$strat{$class};
	$currep = $rep{$class};
	print "$curstrat... $currep...$class \n";
      } else {			#some sort of local use class
	# which is all we really expect to work with templates right now anyway
	$isconstructor=1;
	$issim=$ismessage=0;
	push(@body,$thisline);
      }
    } else {			#classes otherwise unknown
      $issim=$ismessage=0;
      push(@body,$thisline);
    }
    #	print "sim $issim message $ismessage constructor $isconstructor \n";
  } elsif(!$inbody){
    #stuff otherwise unknown
    $outhhandle->print("$thisline\n");
#    $outhhandle->print($lastline{$inhhandle}."\n");
  }
  else
    {
#    print("inbody is $inbody, and we have $thisline\n");
    }
}

$outhhandle->print("#endif\n");
$outhhandle->close;
$inhhandle->close;
# END WRITE .H FILES

# BEGIN WRITE .C FILES
print "Generating source...\n";
my $outChandle=new FileHandle();
if ( $wantindent && ( $indent =~ /GNU/) && ( $indent !~ /2.2.5/) ) {
  $outChandle->open("|indent -npcs -npsl -bap -bad -br -nce>$outC") or die "cannot open $outC";
} else {
  $outChandle->open(">$outC") or die "cannot open $outC";
}

# write outC
$outChandle->print( "#include \"$infile\_sim.h\"\n");
$outChandle->print( "#include \"$infile.def.h\"\n\n");
my $inChandle=new FileHandle();
$inChandle->open("$inC") or die "cannot open $inC";
inithandle($inChandle);
my ($thefnidx);
my($iseventmethod,$isweird,$isnoneventmethod);
my ($returntype,,$messagename,$found);
while (@line=split(' ',($thisline=getcodeline($inChandle)))) {
  if (($line[0] eq "#define") || ($line[0] eq "#include")) {
    $outChandle->print("$thisline\n");
  } elsif ($inbody) {
    if (posefunctest($thisline)>0) {
      push(@body,posefuncmap($thisline,$issim,$isconstructor,$class));
    } else {
      push(@body,$lastline{$inChandle});
    }
  }

  if (($blocklevel{$inChandle}<=1) && ($openblock{$inChandle}==1) &&(!$inbody)) { # decl line
    # handle wrapper class (constructor and resolve) construction for
    # sim classes and event methods
    # get class name and method name
#    	print "open on $thisline :level".$blocklevel{$inChandle}." open:".$openblock{$inChandle}." close:".$closeblock{$inChandle}."\n";
    $iseventmethod=$isnoneventmethod=$isconstructor=$isweird=$issim=0;
    $class=$message=$messagename=$returntype=$method='';
    $declaration=$lastline{$inChandle};
    $inbody=1;
    # regexp if block to extract the fields
    # expecting  class::method(parameters)
    if ($thisline =~ /^\s*([^:\s\(]+)\s*\:\:\s*([^\s\(]+)\s*\(\s*([^\s]+)\s?[\*]*\s*([^\s\,\{\)]+)/) {
      $returntype='';
      $class=$1;
      $method=$2;
      $message=$3;
      $messagename=$4;
    }
    # expecting retval *class::method(parameters)
    elsif ($thisline =~ /^\s*([^:\s\(]+ \*)([^\s:\(]+)\s*\:\:\s*([^\s\(]+)\s*\(\s*([^\s]+)\s?\**\s*([^\s,\)\{]+)/) {

      $returntype=$1;
      $class=$2;
      $method=$3;
      $message=$4;
      $messagename=$5;
    }
    # expecting retval class::method(parameters)
    elsif ($thisline =~ /^\s*([^:\s\(]+)\s?([^\s:\(]+)\s*\:\:\s*([^\s\(]+)\s*\(\s*([^\s]+)\s?\**\s*([^\s,\)\{]+)/) {
      $returntype=$1;
      $class=$2;
      $method=$3;
      $message=$4;
      $messagename=$5;
    }

    # expecting  class::method()
    elsif ($thisline =~ /^\s*\**\s*([^\s:\(]+)\s*::\s*([^\s\(]+)\s*\(\s*\)/) {
      $returntype='';
      $class=$1;
      $method=$2;
      $message=$messagename='';
    }
    # expecting  retval *class::method()
    elsif ($thisline =~  /^\s*([^\s\(:]+\s*\*)([^\s:]+)\s*::\s*([^\s\(]+)\s*\(\s*\)/) {
      $returntype=$1;
      $class=$2;
      $method=$3;
      $message=$messagename='';
    }
    # expecting  retval class::method()
    elsif ($thisline =~  /^\s*([^\s\(:]+)\s?([^\s:]+)\s*::\s*([^\s\(]+)\s*\(\s*\)/) {
      $returntype=$1;
      $class=$2;
      $method=$3;
      $message=$messagename='';
    } else {
      #funkulation: non class subroutine perhaps just pass
      #these for now since we have no handling instructions.
      print STDERR "$thisline at $inC ".$inChandle->input_line_number." is weird \n";
      $isweird=1;
      $inbody=1;
      $outChandle->print($lastline{$inChandle}."\n");
      $declaration='';
      next;
    }
    $issim=1 if exists($base{$class});
    $isconstructor =1 if($class eq $method);
    my $j;
    foreach $j (@{$events{$class}}) {
      if ($method eq $j->[0]) {
	$iseventmethod = 1;
      }
    }
    foreach $j (@{$nonevents{$class}}) {
      if ($method eq $j) {
	$isnoneventmethod = 1;
      }
    }
    #      print "class is $class, method is $method close=".$closeblock{$inChandle}." issim = $issim iseventmethod = $iseventmethod isnoneventmethod $isnoneventmethod \n";

    if ($iseventmethod &&  ($messagename ne '')) {

      $outChandle->print("$declaration\n");
      $outChandle->print("#ifdef POSE_STATS_ON\n");
      $outChandle->print("  int tstat = localStats->TimerRunning();\n");
      $outChandle->print("  if (tstat)\n");
      $outChandle->print("    localStats->SwitchTimer(SIM_TIMER);\n");
      $outChandle->print("  else\n");
      $outChandle->print("    localStats->TimerStart(SIM_TIMER);\n");
      $outChandle->print("#endif\n");
      $outChandle->print("  PVT *pvt = (PVT *)CkLocalBranch(ThePVT);\n");
      $outChandle->print("  Event *e = new Event();\n");
      $outChandle->print("  if ((POSE_endtime < 0) || ($messagename->timestamp <= POSE_endtime)) {\n");
      $outChandle->print("    e->evID = $messagename->evID;\n");
      $outChandle->print("    e->timestamp = $messagename->timestamp;\n");
      $outChandle->print("    e->done = e->commitBfrLen = 0;\n");
      $outChandle->print("    e->commitBfr = NULL;\n");
      $outChandle->print("    e->msg = $messagename;\n");
      $thefnidx = 0;
      #setup index for non constructor methods
      foreach my $i (@{$methods{$class}}) {
	if ($i->[0] ne $class) {
	  foreach my $j (@{$events{$class}}) {

	    if (($i->[0] eq $j->[0])&& ($i->[1] eq $j->[1])) {
	      $thefnidx++;
	    }
	    else
	      {

	      }
	  }
	}
	if (($i->[0] eq $method) && ($i->[1] eq $message)) {
	  print "$class $method $message has index $thefnidx\n";
	  last;
	}
      }
      $outChandle->print("    e->fnIdx = $thefnidx;\n");
      $outChandle->print("    e->next = e->prev = NULL;\n");
      $outChandle->print("    e->spawnedList = NULL;\n");
      $outChandle->print("    int rb = eq->InsertEvent(e);\n");
      $outChandle->print("    if ((e->timestamp < eq->currentPtr->timestamp)\n");
      $outChandle->print("        || (eq->currentPtr == eq->backPtr)) {\n");
      $outChandle->print("      myStrat->ResetRBevent(e);\n");
      $outChandle->print("    }\n");
#      $outChandle->print("    Step();\n");
      $outChandle->print("    if (!rb) Step();\n");
      $outChandle->print("    else if (active == 0) { active = 1; POSE_Objects[thisIndex].Step(); }\n");
      $outChandle->print("  }\n");
      $outChandle->print("  pvt->objUpdate($messagename->timestamp, RECV);\n");
      $outChandle->print("  srVector[$messagename->fromPE]++;\n");
      $outChandle->print("#ifdef POSE_STATS_ON\n");
      $outChandle->print("  if (tstat)\n");
      $outChandle->print("    localStats->SwitchTimer(tstat);\n");
      $outChandle->print("  else\n");
      $outChandle->print("    localStats->TimerStop();\n");
      $outChandle->print("#endif\n");
      $outChandle->print("}\n");
    } elsif ($isnoneventmethod) {
      my $retval =$returntype;
      $retval =undef if $retval =~ /void/;
      $outChandle->print("$declaration\n");
      $outChandle->print("  ".$retval." result;\n ") if($retval);
      $outChandle->print("  result= ") if($retval);
      $outChandle->print("  ((state_".$class." *)objID)->".$method."(".$messagename.");\n");
      $outChandle->print("  return result\n;") if($retval);
      $outChandle->print("}\n");
    } elsif ($issim && $isconstructor &&($messagename ne '')) { 
      #create the wrapper parent class constructor
      $outChandle->print($returntype." ") if($returntype);
      $outChandle->print(join('',$class,"::",$method,'(',$message,' *',$messagename,"){\n"));
      $outChandle->print("#ifdef POSE_STATS_ON\n  localStats->TimerStart(SIM_TIMER);\n#endif\n");
      $outChandle->print("  PVT *pvt = (PVT *)CkLocalBranch(ThePVT);\n");
      $outChandle->print("#ifdef LB_ON  \n");
      $outChandle->print("  LBgroup *localLBG = TheLBG.ckLocalBranch();\n\n");
      $outChandle->print("#endif  \n");
      $outChandle->print("  myStrat = new ".$strat{$class}."();\n");
      $outChandle->print("  $messagename->parent = this;\n");
      $outChandle->print("  $messagename->str = myStrat;\n");
      $outChandle->print("  int _ts = $messagename->timestamp;\n");
      $outChandle->print("#ifdef POSE_STATS_ON\n  localStats->SwitchTimer(DO_TIMER);\n#endif\n");
      $outChandle->print("  objID = new state_$method($messagename);\n");
      $outChandle->print("#ifdef POSE_STATS_ON\n  localStats->SwitchTimer(SIM_TIMER);\n#endif\n");
      $outChandle->print("  myStrat->init(eq, objID, this, thisIndex);\n");
      $outChandle->print("#ifdef POSE_STATS_ON\n  localStats->TimerStop();\n#endif\n");
      $outChandle->print("  myPVTidx = pvt->objRegister(thisIndex, _ts, sync, this);\n");
      $outChandle->print("#ifdef LB_ON  \n");
      $outChandle->print("  myLBidx = localLBG->objRegister(thisIndex, sync, this);\n");
      $outChandle->print("#endif  \n");
      $outChandle->print("}\n");
      #create the pup constructor
      if (exists($puppy{$class})) {
	$outChandle->print(join('','void ',$class,"::pup(PUP::er &p)\n"));
	$outChandle->print("  {\n");
	$outChandle->print("    sim::pup(p);\n");
	$outChandle->print("    if (p.isUnpacking()) {\n");
	$outChandle->print("      myStrat = new ".$strat{$class}.";\n");
	$outChandle->print("      objID = new state_$class(this, myStrat);\n");
	$outChandle->print("      myStrat->init(eq, objID, this, thisIndex);\n");
	$outChandle->print("    }\n");
	$outChandle->print("    ((state_$class *)objID)->pup(p);\n");
	if ($rep{$class} eq "chpt") {
	  $outChandle->print("    Event *ev = eq->front()->next;\n");
	  $outChandle->print("    int checkpointed;\n\n");
	  $outChandle->print("    while (ev != eq->back()) {\n");
	  $outChandle->print("      if (p.isUnpacking()) {\n");
	  $outChandle->print("        p(checkpointed);\n");
	  $outChandle->print("        if (checkpointed) {\n");
	  $outChandle->print("          ev->cpData = new state_$class(this, myStrat);\n");
	  $outChandle->print("          ((state_$class *)ev->cpData)->cpPup(p);\n");
	  $outChandle->print("        }\n");
	  $outChandle->print("        else ev->cpData = NULL;\n");
	  $outChandle->print("      }\n");
	  $outChandle->print("      else {\n");
	  $outChandle->print("        if (ev->cpData) {\n");
	  $outChandle->print("          checkpointed = 1;\n");
	  $outChandle->print("          p(checkpointed);\n");
	  $outChandle->print("          ((state_$class *)ev->cpData)->cpPup(p);\n");
	  $outChandle->print("        }\n");
	  $outChandle->print("        else {\n");
	  $outChandle->print("          checkpointed = 0;\n");
	  $outChandle->print("          p(checkpointed);\n");
	  $outChandle->print("        }\n");
	  $outChandle->print("      } \n");
	  $outChandle->print("     ev=ev->next; \n");
	  $outChandle->print("    }\n");

	}
	$outChandle->print("  }\n");
      }
    } else {			#non constructor
    }
    if ($closeblock{$inChandle}==1) { #its a one liner
      print "oneliner for  $class $method sim $issim\n";
      $declaration =~ s/(\b$class\b)/state_$1/gm      if($issim);
      if (($method eq 'pup') &&(exists($puppy{$class})) && !exists($cppuppy{$class}) && ($rep{$class} eq "chpt")) {
	my $cpdec=$declaration;
	$cpdec =~ s/:pup/:cpPup/gm;
	$outChandle->print($cpdec."\n");
      }
      $outChandle->print($declaration."\n");
      $inbody=0;
    }
    #	print "done open on $thisline :level".$blocklevel{$inChandle}." open:".$openblock{$inChandle}." close:".$closeblock{$inChandle}."\n";
  } elsif (($blocklevel{$inChandle}==0) && ($closeblock{$inChandle}==1)) { #got the whole thing
    #handle verbatim and translation copying of block
    #if rep is chpt and the $puppy{$class} hash has no elements make the cpPup copy
    #      print "on close class is $class, method is $method issim = $issim iseventmethod = $iseventmethod \n";
    #	print "close  for $class $method \n";
    $declaration =~ s/(\b$class\b)/state_$1/gm      if($issim);

    if (($method eq 'pup') &&(exists($puppy{$class}))  && !exists($cppuppy{$class}) && ($rep{$class} eq "chpt")) {
      my $cpdec=$declaration;
      $cpdec =~ s/:pup/:cpPup/gm;
      $outChandle->print($cpdec."\n");
      foreach (@body) {
	if ($_ =~ /(chpt<state[^>]+>)/) { #translate to rep
	  $outChandle->print($`."rep".$'."\n");
	} else {
	  $outChandle->print($_."\n");
	}
      }
    }
    $outChandle->print($declaration."\n");
    $outChandle->print("init($messagename);\n") if(length($messagename) && ($issim && $isconstructor));
    while ($_ =shift @body) {
      s/\bthishandle\b/thisIndex/gm;
      $outChandle->print($_."\n");
    }
    $inbody=0;
  } elsif (!$inbody) {		#regular line
#    print "regular line $thisline at $inC isweird \n";
    $outChandle->print($lastline{$inChandle}."\n");
  } else {
  }
}
$outChandle->print("\n");
$inChandle->close;

#generate ResolveFns 
my ($key,$count,$first,@array);
foreach $key (keys %methods) {
  print "key is $key\n";
  if ($strat{$key} ne "none") {
    my $count = 1;
    my $first = 1;
    my $ifopen=0;
    @array = @{$methods{$key}};
    foreach my $i (@array) {
      $found = 0;
      foreach my $j (@{$events{$key}}) {
	if (($j->[0] eq $i->[0])&& ($i->[1]==$j->[1])){
	  $found = 1;
	}
      }
      if ($found && ($i->[0] ne $key)) {
	if ($ifopen ==0) {
	  $outChandle->print("void $key\:\:ResolveFn(int fnIdx, void *msg)\n{\n");
	  $outChandle->print("  if (fnIdx >0){\n");
	  $outChandle->print("    if (sync == OPTIMISTIC)\n");
	  $outChandle->print("      ((state_$key *) objID)->checkpoint((state_$key *) objID);\n");
	  $outChandle->print("    ((state_$key *) objID)->update(((eventMsg *)msg)->timestamp, ((eventMsg *)msg)->rst);\n");
	  $outChandle->print("  }\n");


	  $outChandle->print("  if (fnIdx == ");
	  $ifopen=1;
	} elsif ($first == 0) {
	  $outChandle->print("  else if (fnIdx == ");
	}
	$outChandle->print("$count) {\n");
	$first = 0;
$outChandle->print("#ifdef POSE_STATS_ON\n");
$outChandle->print("    if (!CpvAccess(stateRecovery)) {localStats->Do();\n");
$outChandle->print("#ifdef POSE_DOP_ON\n");
$outChandle->print("    st = CmiWallTimer();\n");
$outChandle->print("#endif\n");
$outChandle->print("    localStats->SwitchTimer(DO_TIMER);}\n");
$outChandle->print("#endif\n");

	$outChandle->print("    ((state_$key *) objID)->$i->[0](($i->[1] *)msg);\n");
$outChandle->print("#ifdef POSE_STATS_ON\n");
$outChandle->print("    if (!CpvAccess(stateRecovery)) {\n");
$outChandle->print("#ifdef POSE_DOP_ON\n");
$outChandle->print("    et = CmiWallTimer();\n");
$outChandle->print("    eq->currentPtr->ert = eq->currentPtr->srt + (et-st);\n");
$outChandle->print("    ((state_$key *) objID)->ort = eq->currentPtr->ert+0.000001;\n");
$outChandle->print("    eq->currentPtr->evt = ((state_$key *) objID)->OVT();\n");
$outChandle->print("#endif\n");
$outChandle->print("    localStats->SwitchTimer(SIM_TIMER);}\n");
$outChandle->print("#endif\n");

	$outChandle->print("  }\n");
	$outChandle->print("  else if (fnIdx == -$count) {\n");
	$outChandle->print("    ((state_$key *) objID)->$i->[0]_anti(($i->[1] *)msg);\n");
	$outChandle->print("  }\n");
	$count++;
      }
    }
    if ($ifopen) {
      $outChandle->print("}\n\n");
    }
  }
}

foreach $key (keys %methods) {
  if ($strat{$key} ne "none") {
    my $count = 1;
    my $first = 1;
    my $ifopen =0;
    @array = @{$methods{$key}};
    foreach my $i (@array) {
      $found = 0;
      foreach my $j (@{$events{$key}}) {
	if (($j->[0] eq $i->[0]) &&($i->[1]==$j->[1])) {
	  $found = 1;
	}
      }
      if ($found && ($i->[0] ne $key)) {
	if ($ifopen == 0) {
	  $outChandle->print("void $key\:\:ResolveCommitFn(int fnIdx, void *msg)\n{\n");

	  $outChandle->print("  if (fnIdx == ");
	  $ifopen=1;
	} elsif ($first == 0) {
	  $outChandle->print("  else if (fnIdx == ");
	}
	$outChandle->print("$count) {\n");
	$first = 0;
	$outChandle->print("    ((state_$key *) objID)->$i->[0]_commit(($i->[1] *)msg);\n");
	$outChandle->print("  }\n");
	$count++;
      }
    }
    if ($ifopen) {
      $outChandle->print("}\n\n");
    }
  }
}

#create the inline mapsize function
#$outChandle->print("int MapSizeToIdx(int size)\n{\n");
#my $sizeline=0;
#$outChandle->print("    static int eventMsgSz = sizeof(eventMsg);\n"); 
#foreach my $i (keys %emessages) {
#  $outChandle->print("    static int ".$i."Sz = sizeof(".$i.");\n"); 
#}
#$outChandle->print("\n");
#$outChandle->print("    if (size == eventMsgSz) return ".$sizeline++.";\n");
#foreach my $i (keys %emessages) {
#  $outChandle->print("    else if (size == ".$i."Sz) return ".$sizeline++.";\n");
#}
#$outChandle->print("    return $sizeline;\n");
#$outChandle->print("}\n");
#$outChandle->close;

# END WRITE .C FILES


#codelines have  ; or { or */
#we keep reading till our line meets those criteria.
#it could conceivably exceed them and have more than one logical code line
#we're willing to live with that for now. 
#we'll record block open, close, and comment for test convenience.
sub inithandle
  {
    my($inhandle)=@_;
    $lastline{$inhandle}='';
    $openblock{$inhandle}=0;
    $closeblock{$inhandle}=0;
    $blocklevel{$inhandle}=0;
    $opencomment{$inhandle}=0;
    $incomment{$inhandle}=0;
    $closecomment{$inhandle}=0;
  }
#given a file handle read it till we have a proper code line.
#something with either a { } or ; on it.
#Note: we do not handle stuff inside quotes yet.
#throw away all blank lines and comments
#then return the codeline
#keep the original line with comments in lastline 
#so the parser can use the pure uncommented version
#and we can printout the commented version
sub getcodeline
  {
    my ($inhandle)=@_;
    my $retline='';
    my $thisline;
    $lastline{$inhandle}='';
    while ($thisline=$inhandle->getline()) {
      chomp $thisline;
      if ($thisline =~ /^\s*$/m) { #ignore blank lines
	$thisline='';
	next;
      }
      if ($thisline =~ m://.*$:) { # strip single line comment
	$lastline{$inhandle}.=$thisline."\n";
	$thisline=$PREMATCH.$POSTMATCH;
      } else {
	$lastline{$inhandle}.=$thisline;
      }
      $retline.=$thisline;

      if ($retline =~ m:^\s*/\*.*\*/\s*$:m ) {
	#pure comment line
	# throw it out
	$retline='';
	$closecomment{$inhandle}=0;
	$opencomment{$inhandle}=0;
	$incomment{$inhandle}=0;
	$lastline{$inhandle}.="\n";
	next;
      } elsif ($retline =~ m:/\*.*\*/:m) { # line containing comment plus other stuff
	#strip comments out
	while ($retline =~ m:/\*.*\*/:m) {
	  $retline=$PREMATCH.$POSTMATCH;
	  $closecomment{$inhandle}=0;
	  $incomment{$inhandle}=0;
	  $opencomment{$inhandle}=0;
	}
	$lastline{$inhandle}.="\n";
	next if($retline =~ /^\s*$/); #nothing to test
      } elsif ($retline =~ m:/\*:m) {
	$lastline{$inhandle}.="\n";
	# we don't really use opencomment or closecomment for
	# anything, but might as well keep track of the state.
	$closecomment{$inhandle}=1;
	next; #should now be handled by the comment stripper condition
      } elsif ($retline =~ m:\*/:m) {
	$opencomment{$inhandle}=1;
	$incomment{$inhandle}=1;
	next;	     #keep scanning till we got the end of the comment
      } elsif ($incomment{$inhandle}==1) { #keep scanning for end
	next;
      }
      if ($retline =~ /\{/m) {
	$openblock{$inhandle}=1;
	$blocklevel{$inhandle}+=($retline =~tr/\{//);
      } else {
	$openblock{$inhandle}=0;
      }
      if ($retline =~ /\}/m) {
	$closeblock{$inhandle}=1;
	$blocklevel{$inhandle}-=($retline =~tr/\}//);
      } else {
	$closeblock{$inhandle}=0;
      }
      last if ($retline =~ /^\#/);
      last if ($retline =~ /;/m);
      last if ($retline =~ /^\s*\w+\s*:\s*$/); # let publi|cprivate on its own line pass through 
      last if($openblock{$inhandle}==1);
      last if($closeblock{$inhandle}==1);
    }
    return $retline;
  }

#given first line of a char object
#record useful structure
sub poserdeal
  {
    my($inhandle,$lineref,$thisline)=@_;
    my $poser;
    my @line=@$lineref;
    my $class;
    if (($openblock{$inhandle}==1) && ($blocklevel{$inhandle}==1)) {
      if ($thisline =~ /^\s*(poser)\s+(\S+)\s+\:\s+(sim)\s+(\S+)\s+(\S+)\s*\{\s*$/) {

	$line[0]="chare";
	$class=$2;
	print " poserdeal class $class \n";
	push(@strategies, $4);
	push(@representations , $5);
	$base{$class}=$3;
	$strat{$class}=$4;
	$rep{$class}=$5;
	$events{$class}=[];
	push(@{$poser->{outarr}}, "  @line\n");
	push(@{$poser->{literal}}, "array [1D] ".$class." : ".$3." \{");
	print  "chare ".$class." : ".$3."\{\n";
	$needsopeq{$class}=1 if($strat{$class} eq "opt");
	$needsopeq{$class}=1 if($rep{$class} eq "chpt");
      } else {
	print $thisline;
	die "bad poser declaration";
      }
    } else {
      print $thisline;
      die "bad poser declaration";
    }
    #if sim, add to the base strat and rep hashes
    #scan for event declarations, track them in events hash and strip event keyword
    #store non even entries
    while (@line=split(' ',($thisline=getcodeline($inhandle)))) {
      if ($closeblock{$inhandle}==1) {
	push(@{$poser->{literal}}, $lastline{$inhandle});
	push(@{$poser->{outarr}}, "  ".$thisline);
	last;
      } elsif ($line[0] eq "entry") {
	if ($thisline =~ /(\,event\,|\,event|event\,|\s\[event\]).*\s(.*)\(\s*(\S[^\*\s\)]*)/) {
	  $method = $2;
	  $message = $3;
	  push(@{$events{$class}}, [$method,$message]);
#	  print " poserdeal events  class $class  method $method message $message \n";
	} elsif ($thisline =~ /\s+[\*]*\s*(\S+)\(/) {
	  $method = $1;
	  push(@{$nonevents{$class}}, ($method)) if $class ne $method;
	}
	$thisline =~ s/(\,event\,|\,event|event\,|\s\[event\])//;
	push(@{$poser->{literal}}, $thisline);
	push(@{$poser->{entry}},$thisline);
	push(@{$poser->{outarr}}, "  ".$thisline);
      }
    }
    return $poser;
  }


sub groupnchardeal
  {
    my($inhandle,$lineref,$thisline)=@_;
    my $group;
    my @line=@$lineref;
    my $class;
    push(@{$group->{literal}}, $lastline{$inhandle});
    if (($line[2] eq "{") || ($line[4] eq "{")) {
      $strat{$line[1]} = "none";
      push(@{$group->{outarr}}, "  $thisline\n");

    } else {
      die "what the hell is this? $thisline ";
    }
    #if sim, add to the base strat and rep hashes
    #scan for event declarations, track them in events hash and strip event keyword
    while ($thisline=getcodeline($inhandle)) {
      @line=split(' ',$thisline);
      push(@{$group->{literal}}, $lastline{$inhandle});
      if ($closeblock{$inhandle}) {
	push(@{$group->{outarr}}, "  ".$thisline);
	last;
      } elsif ($line[0] eq "entry") {
	if ($thisline =~ /(\,event\,|\,event|event\,|\s\[event\]).*\s(.*)\(\s*(\S[^\*\s\)]*)/) {
#	  $method = $2;
	  $message = $3;
	  $method = $line[2];
	  $method =~ s/\*?(.+)\(.+/$1/;
	  push(@{$events{$class}}, [$method,$message]);
	} elsif ($thisline =~ /\s+[\*]*\s*(\S+)\(/) {	
	  $method = $1;
	  push(@{$nonevents{$class}}, ($method)) if $class ne $method;
	}
	$thisline =~ s/(\,event\,|\,event|event\,|\s\[event\])//;
	push(@{$group->{entry}},$thisline);
	push(@{$group->{outarr}}, "  ".$thisline);
      } else {
	if (($line[2] eq "{") || ($line[4] eq "{")) {
	  $strat{$line[1]} = "none";
	  push(@{$group->{outarr}}, "  ".$thisline);
	} else {
	  die "what the hell is this $thisline";
	}
      }
    }
    return $group;
  }

sub posefunctest
  {
    my ($line)=shift;
    if ($line =~ /POSE_create/) {
      return 1;
    } elsif ($line =~ /POSE_creations_complete/) {
      return 2;
    } elsif ($line =~ /POSE_invoke_at/) {
      return 4;
    } elsif ($line =~ /POSE_invoke/) {
      return 3;
    } elsif ($line =~ /POSE_local_invoke/) {
      return 5;
    } else {
      return 0;
    }
  }

#given a POSE_* call break it down into its segments and return them as a list
sub posefuncparse
  {
    my($line)=@_;
    my($restline);
    my @segments;
    if ($line =~ /\)\s*;\s*$/) { #strip );
      my $foo=$MATCH;
      $foo =~ tr/ \t//d;
      $line=$PREMATCH.$foo;
      chop $line;
      chop $line;
    }
    if ($line =~ /(POSE_[\w_]+)\s*\(/) { #strip funcname(
      $segments[0]=$1;
    } else {
      return undef;
    }
    my $preline=$PREMATCH;
    $restline=$POSTMATCH;
    if (($restline =~/\(+[^\)]*,[^\)]*\)/) || ($restline =~/\[+[^\]]*,[^\]]*\]/)) {
      #ok now we hunt down the groups and take them out of the expression
      #replace with &*# number
      my @inner_ugliness;
      my $hackline=$restline;
      while (($hackline =~/(\([^\(\)]*,[^\(\)]*\))/)||($hackline =~/(\[[^\[\]]*,[^\]\[]*\])/)) { # there exists an inner parenthesized expression with commas in
	push(@inner_ugliness,$1);
	$hackline=$PREMATCH.' &*#'.$#inner_ugliness.' '.$POSTMATCH;
      }
      my @uglysegs=split(/,/,$hackline);
      my $i;
      for ($i=0,$i<=$#uglysegs,$i++) {
	while ($uglysegs[$i] =~ /( &\*\#)([\d]+ )/) {
	  $uglysegs[$i]=$PREMATCH.$inner_ugliness[$2].$POSTMATCH;
	}
      }
      push(@segments,@uglysegs);
    } else {			#no ickiness just split it
      push(@segments,split(/,/,$restline));
    }
    my $i;
    for ($i=0;$i<$#segments;$i++) {
      $segments[$i]=~ s/^\s*(.*)\s*$/$1/e;
    }
    return($preline,@segments);
  }

#given an input line and the context
#determine if it is a posefunc and if so map it.
sub posefuncmap
  {
    my($line,$issim,$isconstructor,$simobjtype)=@_;
    my($type,$output);
    if (($type=posefunctest($line))>0) {
      #do stuff
      my($preline,@segments)=posefuncparse($line);
      if ($type==1)		#create
	{
	  my($sim,$msg);
	  if ($segments[1]=~/\s*([^(]*)\s*\(\s*([^)]*)\s*\)/) {
	    $sim=$1;
	    $msg=$2;
	    $output=$preline."\n{\n";
	    $output.="int _POSE_handle = ".$segments[3].";\n";
	    $output.="int _POSE_atTime = ".$segments[4].";\n" if ($#segments>=4);
	    $output.=$msg."->Timestamp(_POSE_handle);\n";
	    $output.="$msg->rst = 0.0;\n";

	    $output.="(*(CProxy_".$sim." *)&POSE_Objects)[".$segments[2]."].insert(".$msg;
	    if ($#segments>=4) {
	      $output.=",_POSE_atTime";
	    }
	    $output.=");}\n";
	  } else {
	    die "what the hell is up with .".$segments[1]."in $line";
	  }
	} elsif ($type==2)	#create_complete
	  {
	    $output=$preline."POSE_Objects.doneInserting();\n";
	  } elsif ($type==3)	#invoke
	    {
	      if (!$issim) {
		die "cannot use POSE_invoke from non sim object";
	      }
	      my($event,$msg);
	      if ($segments[1]=~/\s*([^(]*)\s*\(\s*([^)]*)\s*\)/) {
		$event=$1;
		$msg=$2;
		if ($isconstructor) {
		  $output=$preline."\n{\n";
		  $output.="int _POSE_handle = ".$segments[3].";\n";
		  $output.="unsigned int _POSE_timeOffset = ".$segments[4].";\n" if ($#segments>=4);
		  $output.="PVT *pvt = (PVT *)CkLocalBranch(ThePVT);\n";
		  $output.=$msg."->Timestamp(ovt+(_POSE_timeOffset));\n";
		  $output.="pvt->objUpdate(ovt+(_POSE_timeOffset), SEND);\n";
		  $output.="$msg->fromPE = CkMyPe();\n";
		  $output.="$msg->rst = 0.0;\n";
		  $output.="(* (CProxy_".$segments[2]." *)&POSE_Objects)[_POSE_handle].".$segments[1].";\n";
		  $output.="int _destPE = POSE_Objects.ckLocalBranch()->lastKnown(CkArrayIndex1D(_POSE_handle));\n";
		  $output.="parent->srVector[_destPE]++;\n";
		  $output.="}\n";
		} else {
		  $output=$preline."\n{\n";
		  $output.="if (!CpvAccess(stateRecovery)) {\n";
		  $output.="int _POSE_handle = ".$segments[3].";\n";
		  $output.="unsigned int _POSE_timeOffset = ".$segments[4].";\n";
		  $output.="registerTimestamp(_POSE_handle, ".$msg.",_POSE_timeOffset);\n";
		  $output.="$msg->fromPE = CkMyPe();\n";
		  $output.="#ifdef POSE_DOP_ON\n";
		  $output.="parent->ct = CmiWallTimer();\n";
		  $output.="$msg->rst = parent->ct - parent->st + parent->eq->currentPtr->srt;\n";
		  $output.="#endif\n";
		  $output.="(* (CProxy_".$segments[2]." *)&POSE_Objects)[_POSE_handle].".$segments[1].";\n";
		  $output.="int _destPE = POSE_Objects.ckLocalBranch()->lastKnown(CkArrayIndex1D(_POSE_handle));\n";
		  $output.="parent->srVector[_destPE]++;\n";
		  $output.="}\n";
		  $output.="else delete ".$msg.";}\n";
		}
	      } else {
		die "what the hell is up with .".$segments[1]."in $line";
	      }
	    } elsif ($type==4)	#invoke_at
	      {
		my($event,$msg);
		if ($segments[1]=~/\s*([^(]*)\s*\(\s*([^)]*)\s*\)/) {
		  $event=$1;
		  $msg=$2;
		  if (!$issim || ($issim && $isconstructor)) {
		    print "warning: should use POSE_invoke instead of POSE_invoke_at inside sim object constructors\n" if ($issim && $isconstructor);
		    $output=$preline."\n{\n";
		    $output.="int _POSE_handle = ".$segments[3].";\n";
		    $output.="int _POSE_atTime = ".$segments[4].";\n";
		    $output.="PVT *pvt = (PVT *)CkLocalBranch(ThePVT);\n";
		    $output.=$msg."->Timestamp(_POSE_atTime);\n";
		    $output.="pvt->objUpdate(_POSE_atTime, SEND);\n";
		    $output.="$msg->fromPE = CkMyPe();\n";
		  $output.="$msg->rst = 0.0;\n";
		    $output.="(*(CProxy_".$segments[2]." *)&POSE_Objects)[_POSE_handle].".$segments[1].";\n";
		    $output.="int _destPE = POSE_Objects.ckLocalBranch()->lastKnown(CkArrayIndex1D(_POSE_handle));\n";
		    $output.="parent->srVector[_destPE]++;\n";
		    $output.="}\n";
	      
		    $output.="}\n";
		  } else {
		    print "warning: should use POSE_invoke instead of POSE_invoke_at in side sim object nonconstructors\n" if ($issim && $isconstructor);
		    $output=$preline."\n{\n";
		    $output.="int _POSE_handle = ".$segments[3].";\n";
		    $output.="int _POSE_atTime = ".$segments[4].";\n";
		    $output.="if (!CpvAccess(stateRecovery)) {\n";
		    $output.="registerTimestamp(_POSE_handle, ".$msg.", _POSE_atTime);\n";
		    $output.="$msg->fromPE = CkMyPe();\n";
		  $output.="#ifdef POSE_DOP_ON\n";
		  $output.="parent->ct = CmiWallTimer();\n";
		  $output.="$msg->rst = parent->ct - parent->st + parent->eq->currentPtr->srt;\n";
		  $output.="#endif\n";
		    $output.="(* (CProxy_".$segments[2]." *)&POSE_Objects)[_POSE_handle].".$segments[1].";\n";
		    $output.="int _destPE = POSE_Objects.ckLocalBranch()->lastKnown(CkArrayIndex1D(_POSE_handle));\n";
		    $output.="parent->srVector[_destPE]++;\n";
		    $output.="}\n";
		    $output.="else delete ".$msg.";}\n";
		  }
		} else {
		  die "what the hell is up with .".$segments[1]."in $line";
		}
	      } elsif ($type==5) #local_invoke
		{
		  my($event,$msg);
		  if ($segments[1]=~/\s*([^(]*)\s*\(\s*([^)]*)\s*\)/) {
		    $event=$1;
		    $msg=$2;
		    if ($issim && $isconstructor) {
		      $output=$preline."\n{\n";
		      $output.="unsigned int _POSE_timeOffset = ".$segments[2].";\n";
		      $output.="PVT *pvt = (PVT *)CkLocalBranch(ThePVT);\n";
		      $output.=$msg."->Timestamp(ovt+(_POSE_timeOffset));\n";
		      $output.="pvt->objUpdate(ovt+(_POSE_timeOffset), SEND);\n";
		      $output.="$msg->fromPE = CkMyPe();\n";
		  $output.="$msg->rst = 0.0;\n";
		      $output.="(* (CProxy_".$simobjtype." *)&POSE_Objects)[parent->thisIndex].".$segments[1].";\n";
		      $output.="parent->srVector[CkMyPe()]++;\n";
		      $output.="}\n";
		    } elsif ($issim) {
		      $output=$preline."\n{\n";
		      $output.="unsigned int _POSE_timeOffset = ".$segments[2].";\n";
		      $output.="if (!CpvAccess(stateRecovery)) {\n";
		      $output.="registerTimestamp(parent->thisIndex, ".$msg.", _POSE_timeOffset);\n";
		      $output.="$msg->fromPE = CkMyPe();\n";
		  $output.="#ifdef POSE_DOP_ON\n";
		  $output.="parent->ct = CmiWallTimer();\n";
		  $output.="$msg->rst = parent->ct - parent->st + parent->eq->currentPtr->srt;\n";
		  $output.="#endif\n";
		      $output.="(* (CProxy_".$simobjtype." *)&POSE_Objects)[parent->thisIndex].".$segments[1].";\n";
		      $output.="parent->srVector[CkMyPe()]++;\n";
		      $output.="}\n";
		      $output.="else delete(".$msg.");\n}";
		    } else {
		      die "untranslatably deranged to call POSE_local_invoke from non sim object";
		    }
		  } else {
		    die "what the hell is up with .".$segments[1];
		  }
		}
    }
  }

