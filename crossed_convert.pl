#!/usr/bin/perl -w

#CRoSSeD v1.1

use warnings;
use strict;
use POSIX qw(floor ceil);

#fasta-format input file
my$fastafile = "fasta.txt";
if (defined($ARGV[0])){
	$fastafile = $ARGV[0];
}

#crossed/crf-format output file
my$outputfile = "convertedfasta.txt";
if (defined($ARGV[1])){
	$outputfile = $ARGV[1];
}

#define type of output file (1: crossed screening, 2: standard crf++ input)
my$type = 2;
if (defined($ARGV[2])){
	$type = $ARGV[2];
}

my$defaultlabel;
if (defined($ARGV[3])){
	$defaultlabel = $ARGV[3];
}

open(FASTA,$fastafile) or die "Unable to open fasta file\n";
my$nr_seq = -1;
my@data;
while(<FASTA>){
	chomp;
	my$string = $_;
	
	if($string =~ m/^>/){
		if ($nr_seq == -1){
			$nr_seq++;
		} else {
			if(defined($data[$nr_seq])){
				$nr_seq++;
			}
		}
	} elsif ($string =~ m/^(A|T|G|C)/i){
		if (defined $data[$nr_seq]){
			$data[$nr_seq] = $data[$nr_seq].$string;
		} else {
			$data[$nr_seq] = $string;
		}
	}
}
close FASTA;

open (OUTPUT, ">".$outputfile) or die "Cannot open output file\n";
foreach my$sequence (@data){
	$sequence =~ s/(\w)/\U$1/gi;
	$sequence =~ s/\W+//gi;
	my@tempmat;
	push @tempmat, [split("",$sequence)];
	push @tempmat, [aphil(@{$tempmat[0]})];
	push @tempmat, [bend(@{$tempmat[0]})];
	push @tempmat, [btwist(@{$tempmat[0]})];
	push @tempmat, [deform(@{$tempmat[0]})];
	push @tempmat, [denat(@{$tempmat[0]})]; #Deform used to be switched with denat in some versions
	push @tempmat, [disrupt(@{$tempmat[0]})];
	push @tempmat, [gc(@{$tempmat[0]})];
	push @tempmat, [prop(@{$tempmat[0]})];
	push @tempmat, [protwist(@{$tempmat[0]})];
	push @tempmat, [stab(@{$tempmat[0]})];
	push @tempmat, [stack(@{$tempmat[0]})];
	push @tempmat, [stiff(@{$tempmat[0]})];
	push @tempmat, [zdna(@{$tempmat[0]})];
	if ($type == 1){
		if(defined($defaultlabel)){
			print OUTPUT $defaultlabel."\n";
		}
		for my$i (0 .. (length($sequence)-1)){
			print OUTPUT $tempmat[0][$i];
			for my$j (1 .. $#tempmat){
				print OUTPUT "\t".$tempmat[$j][$i];
			}
			print OUTPUT "\n";
		}
		print OUTPUT "\n";
	} elsif ($type == 2) {
		print OUTPUT join("\t",@{$tempmat[0]})."\t";
		for my$j (1 .. $#tempmat){
			print OUTPUT join("\t",@{$tempmat[$j]})."\t";
		}
		if(defined($defaultlabel)){
			print OUTPUT $defaultlabel;
		}
		print OUTPUT "\n";
	}
}
close OUTPUT;



sub aphil {
	my@array = @_;
	
	my%a = ("A",0.97,"T",	0.58,"G",0.33,"C",0.13);
	my%t = ("A",0.73,"T",	0.97,"G",1.04,"C",0.98);
	my%g = ("A",0.98,"T", 0.13,"G",0.19,"C",0.73);
	my%c = ("A",1.04,"T",	0.33,"G",0.52,"C",0.19);
	my%mat =  ("A",\%a,"T",\%t,"G",\%g,"C",\%c);

	my@stfeat;
	for my$x (1 .. $#array){
		push @stfeat, $mat{$array[$x-1]}{$array[$x]};
	}
	unshift @stfeat, '+';

	return @stfeat;
	
}

sub bend {

	my@array = @_;
	my%aa = ("A",  -0.2740,"T",    0.1820,"G" ,   0.0270,"C",   -0.0060);
	my%at = ("A", 0.0680,"T",    0.0680,"G",    0.1940,"C",    0.1940);
	my%ag = ("A", -0.0370,"T",    0.0250,"G",    0.0130,"C",    0.0760);
	my%ac = ("A", 0.0150,"T",    0.0900,"G",   -0.0030,"C",   -0.2460);
	my%a = ("A",\%aa,"T",\%at,"G",\%ag,"C",\%ac);

	my%ta = ("A",    -0.2800,"T",   -0.2800,"G" ,   -0.1830 ,"C",  -0.1830);
	my%tt = ("A",     0.1820 ,"T",  -0.2740 ,"G" ,  -0.0060  ,"C",  0.0270);
	my%tg = ("A",    -0.1100 ,"T",  -0.2050 ,"G" ,  -0.0320 ,"C",   0.0170);
	my%tc = ("A",     0.1340 ,"T",  -0.0810 ,"G" ,  -0.0330 ,"C",  -0.0570);
	my%t = ("A",\%ta,"T",\%tt,"G",\%tg,"C",\%tc);


	my%ga = ("A",  -0.0810,"T",    0.1340,"G" ,   -0.0570 ,"C",  -0.0330);
	my%gt = ("A",   0.0900,"T",    0.0150 ,"G" ,  -0.2460 ,"C",  -0.0030);
	my%gg = ("A",   0.0310 ,"T",   0.0400 ,"G" ,  -0.0120 ,"C",  -0.0770);
	my%gc = ("A",   0.1750,"T",    0.1750 ,"G" ,  -0.1360,"C",  -0.1360);
	my%g = ("A",\%ga,"T",\%gt,"G",\%gg,"C",\%gc);


	my%ca = ("A",   -0.2050,"T",   -0.1100,"G" ,    0.0170,"C",   -0.0320);
	my%ct = ("A",    0.0250 ,"T",  -0.0370,"G" ,    0.0760,"C",    0.0130);
	my%cg = ("A",   -0.0130 ,"T",  -0.0130 ,"G" ,   0.1070,"C",    0.1070);
	my%cc = ("A",    0.0400 ,"T",   0.0310 ,"G" ,  -0.0770 ,"C",  -0.0120);
	my%c = ("A",\%ca,"T",\%ct,"G",\%cg,"C",\%cc);

	my%mat =  ("A",\%a,"T",\%t,"G",\%g,"C",\%c);
	my@stfeat;
	for my$x (2 .. $#array){
		push @stfeat, $mat{$array[$x]}{$array[$x-2]}{$array[$x-1]};
	}
	unshift @stfeat, '+';
	unshift @stfeat, '+';
	
	return @stfeat;
 
}

sub btwist {

	my@array = @_;
	my%a = ("A",35.5,"T",	43.2,"G",	30.6,"C",	33.1);
	my%t = ("A",31.6,"T",	35.5,"G",	37.7,"C",	39.6);
	my%g = ("A",39.6,"T",	33.1,"G",	35.3,"C",	38.4);
	my%c = ("A",37.7,"T",	30.6,"G",	31.3,"C",	35.3);
	my%mat =  ("A",\%a,"T",\%t,"G",\%g,"C",\%c);
	
	my@stfeat;
	for my$x (1 .. $#array){
		push @stfeat, $mat{$array[$x-1]}{$array[$x]};
	}
	unshift @stfeat, '+';
	
	return @stfeat;

}

sub deform {

	my@array = @_;

	my%a = ("A",2.9,"T",	1.6,"G",	2.1,"C",	2.3);
	my%t = ("A",6.3,"T",	2.9,"G",	9.8,"C",	4.5);
	my%g = ("A",4.5,"T",	2.3,"G",	6.1,"C",	4);
	my%c = ("A",9.8,"T",	2.1,"G",	12.1,"C",	6.1);
	my%mat =  ("A",\%a,"T",\%t,"G",\%g,"C",\%c);
	
	my@stfeat;
	for my$x (1 .. $#array){
		push @stfeat, $mat{$array[$x-1]}{$array[$x]};
	}
	unshift @stfeat, '+';
	
	return @stfeat;
}

sub denat {

	my@array = @_;
	my%a = ("A",66.51,"T",	72.29,"G",	85.12,"C",	108.8);
	my%t = ("A",50.11,"T",	66.51,"G",	64.92,"C",	80.03);
	my%g = ("A",80.03,"T",	64.92,"G",	99.31,"C",	135.83);
	my%c = ("A",64.92,"T",	85.15,"G",	88.84,"C",	99.31);
	my%mat =  ("A",\%a,"T",\%t,"G",\%g,"C",\%c);
	
	my@stfeat;
	for my$x (1 .. $#array){
		push @stfeat, $mat{$array[$x-1]}{$array[$x]};
	}
	unshift @stfeat, '+';
	
	return @stfeat;
}

sub disrupt {

	my@array = @_;
	my%a = ("A",1.9,"T",	0.9,"G",	1.6,"C",	1.3);
	my%t = ("A",1.5,"T",	1.9,"G",	1.9,"C",	1.6);
	my%g = ("A",1.6,"T",	1.3,"G",	3.1,"C",	3.1);
	my%c = ("A",1.9,"T",	1.6,"G",	3.6,"C",	3.1);
	my%mat =  ("A",\%a,"T",\%t,"G",\%g,"C",\%c);
	
	my@stfeat;
	for my$x (1 .. $#array){
		push @stfeat, $mat{$array[$x-1]}{$array[$x]};
	}
	unshift @stfeat, '+';
	
	return @stfeat;
}

sub gc {

	my@array = @_;
	my%a = ("A",0,"T",	0,"G",	1,"C",	1);
	my%t = ("A",0,"T",	0,"G",	1,"C",	1);
	my%g = ("A",1,"T",	1,"G",	2,"C",	2);
	my%c = ("A",1,"T",	1,"G",	2,"C",	2);
	my%mat =  ("A",\%a,"T",\%t,"G",\%g,"C",\%c);

	my@stfeat;
	for my$x (1 .. $#array){
		push @stfeat, $mat{$array[$x-1]}{$array[$x]};
	}
	unshift @stfeat, '+';
	
	return @stfeat;

}

sub prop {

	my@array = @_;
	my%a = ("A",-18.66,"T",	-15.01,"G",	-14,"C",	-13.1);
	my%t = ("A",-11.85,"T",	-18.66,"G",	-9.45,"C",	-13.48);
	my%g = ("A",-13.48,"T",	-13.1,"G",	-8.11,"C",	-11.08);
	my%c = ("A",-9.45,"T",	-14	,"G",-10.03,"C",	-8.11);
	my%mat =  ("A",\%a,"T",\%t,"G",\%g,"C",\%c);

	my@stfeat;
	for my$x (1 .. $#array){
		push @stfeat, $mat{$array[$x-1]}{$array[$x]};
	}
	unshift @stfeat, '+';
	
	return @stfeat;
}

sub protwist {

	my@array = @_;
	my%a = ("A",35.1,"T",	29.3,"G",	31.9,"C",	31.5);
	my%t = ("A",37.8,"T",	35.1,"G",	37.3,"C",	36.3);
	my%g = ("A",36.3,"T",	31.5,"G",	32.9,"C",	33.6);
	my%c = ("A",37.3,"T",	31.9,"G",	36.1,"C",	32.9);
	my%mat =  ("A",\%a,"T",\%t,"G",\%g,"C",\%c);

	my@stfeat;
	for my$x (1 .. $#array){
		push @stfeat, $mat{$array[$x-1]}{$array[$x]};
	}
	unshift @stfeat, '+';
	
	return @stfeat;

}

#Sugimoto free energy
sub stab {

	my@array = @_;
	my%a = ("A",-1.2,"T",	-0.9,"G",	-1.5,"C",	-1.5);
	my%t = ("A",-0.9,"T",	-1.2,"G",	-1.7,"C",	-1.5);
	my%g = ("A",-1.5,"T",	-1.5,"G",	-2.1,"C",	-2.3);
	my%c = ("A",-1.7,"T",	-1.5,"G",	-2.8,"C",	-2.1);
	my%mat =  ("A",\%a,"T",\%t,"G",\%g,"C",\%c);

	my@stfeat;
	for my$x (1 .. $#array){
		push @stfeat, $mat{$array[$x-1]}{$array[$x]};
	}
	unshift @stfeat, '+';
	
	return @stfeat;
}


sub stack {

	my@array = @_;
	my%a = ("A",-5.37,"T",	-6.57,"G",	-6.78,"C",	-10.51);
	my%t = ("A",-3.82,"T",	-5.37,"G",	-6.57,"C",	-9.81);
	my%g = ("A",-9.81,"T",	-10.51,"G",	-8.26,"C",	-14.59);
	my%c = ("A",-6.57,"T",	-6.78,"G",	-9.61,"C",	-8.26);
	my%mat =  ("A",\%a,"T",\%t,"G",\%g,"C",\%c);

	my@stfeat;
	for my$x (1 .. $#array){
		push @stfeat, $mat{$array[$x-1]}{$array[$x]};
	}
	unshift @stfeat, '+';
	
	return @stfeat;
}

sub stiff {

	my@array = @_;
	my%a = ("A",35,"T",	20,"G",	60,"C",	60);
	my%t = ("A",20,"T",	35,"G",	60,"C",	60);
	my%g = ("A",60,"T",	60,"G",	130,"C",	85);
	my%c = ("A",60,"T",	60,"G",	85,"C",	130);	
	my%mat =  ("A",\%a,"T",\%t,"G",\%g,"C",\%c);

	my@stfeat;
	for my$x (1 .. $#array){
		push @stfeat, $mat{$array[$x-1]}{$array[$x]};
	}
	unshift @stfeat, '+';
	
	return @stfeat;

}

sub zdna {

	my@array = @_;
	my%a = ("A",3.9,"T",	5.9,"G",	3.4,"C",	4.6);
	my%t = ("A",2.5,"T",	3.9,"G",	1.3,"C",	3.4);
	my%g = ("A",3.4,"T",	4.6,"G",	2.4,"C",	4);
	my%c = ("A",1.3,"T",	3.4,"G",	0.7,"C",	2.4);
	my%mat =  ("A",\%a,"T",\%t,"G",\%g,"C",\%c);

	my@stfeat;
	for my$x (1 .. $#array){
		push @stfeat, $mat{$array[$x-1]}{$array[$x]};
	}
	unshift @stfeat, '+';
	
	return @stfeat;


}