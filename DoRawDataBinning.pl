#!/usr/pubsw/bin/perl

use strict;
use Cwd;
use File::Copy;

my $mainDir = getcwd;
my $binProcFile = "/autofs/space/gulliver_001/users/Baoqiang/capillaryPO2/code4po2/basic PO2/bindata";
my @theFileNames;
my $drc;
my $file;
my $isDir;
my $rawFlag;
my $zippedFlag;
my $processedFlag;
my $status;
my @binout;
my $doneFlag;
my $foo;
my $edgeFlag;
my $eqIncrementFlag;
my $incrWritten;
my $incrCalc;
my @foo1;
my $foo2;
my $idx;
my @tarout;

opendir(IMD,$mainDir) || die("Can't open directory");
my @theDirNames = readdir(IMD);

@theDirNames = sort(@theDirNames);
foreach $drc (@theDirNames){

	if(!(($drc eq ".") || ($drc eq ".."))){
		#print "Oppening $drc \n";
		$isDir = opendir(IMD1,$drc);
		if($isDir == 1){
			#print "Oppened $drc \n";
			@theFileNames = readdir(IMD1);
			$rawFlag = 0;
			$zippedFlag = 0;
			$processedFlag = 0;
			foreach $file (@theFileNames){
				if($file eq "raw_data.bin"){
					$rawFlag = 1;}
				if($file eq "raw_data.bin.tar.gz"){
					$zippedFlag = 1;}
				if($file eq "binned_data.txt"){
					$processedFlag = 1;}
			}
			if($rawFlag && !$zippedFlag){
				# copy program, process, wait, check, delete...
				chdir($drc);
					copy($binProcFile,"bindata");
					chmod(0774,"bindata");
					print "Working in directory $drc \n";
					system("./bindata 2 2 1>bindataout.txt 2>bindataouterr.txt");
#					@binout = `bindata 2 2`;
					open(FID,"bindataout.txt");
					@binout = <FID>;
					close(FID);
#					open FID, ">", "bindataout.txt";
					$doneFlag = 0;
					$edgeFlag = 1;
					$eqIncrementFlag = 0;
					foreach $foo (@binout){
						if($foo =~ m/done./){
							$doneFlag = 1;}
						if($foo =~ m/Edge/){
							$edgeFlag = 0;}
						if($foo =~ m/Binned data file written./){
							@foo1 = split(/ /,$foo);
							foreach $foo2 (@foo1){
								$idx = $idx + 1;
								if($idx eq 5){
									$incrWritten = $foo2;
									#print FID ("Written  $incrWritten \n"); 
								}
								else{
									if($idx eq 8){
										$incrCalc = $foo2;
										#print FID ("Calcul  $incrCalc \n");
									}
								}	
							}
						}
#						print FID $foo;
					}
#					close(FID);
					if($doneFlag && $edgeFlag && ($incrWritten eq $incrCalc)){
						open FID, ">", "success_bin_data";
						close(FID);
						#do tar gz of raw_data.bin
						#@tarout = `tar -cvzf raw_data.bin.tar.gz raw_data.bin`;
					} else{
						open FID, ">", "error_bin_data";
						close(FID);
					}
					
				chdir($mainDir);
			}
			closedir(IMD1);
		}
	}
}

closedir(IMD);


