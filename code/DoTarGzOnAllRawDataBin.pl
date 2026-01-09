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
			if($rawFlag && !$zippedFlag && $processedFlag){
				# copy program, process, wait, check, delete...
				chdir($drc);
				print "Working in directory $drc \n";
				system("tar -cvzf raw_data.bin.tar.gz raw_data.bin 1>tarout.txt 2>tarouterr.txt");
				chdir($mainDir);
			}
			closedir(IMD1);
		}
	}
}

closedir(IMD);


