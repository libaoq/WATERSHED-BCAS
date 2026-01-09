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

print "WARNING!!! This will delete all raw_data.bin files!!! \n";
print "Make sure that you have a backup!!! \n";
print "This function will not delete files unless there is tar.gz archive! \n\n";
print "Are you sure that you want to proceed? (YeS/no):  ";
my $userInpt = <>;
chomp($userInpt);

if ($userInpt eq "YeS"){
	print "\nO.K. deleting started... \n";

	opendir(IMD,$mainDir) || die("Can't open directory");
	my @theDirNames = readdir(IMD);
	@theDirNames = sort(@theDirNames);
	
	foreach $drc (@theDirNames){
		if(!(($drc eq ".") || ($drc eq ".."))){
			$isDir = opendir(IMD1,$drc);
			if($isDir == 1){
				@theFileNames = readdir(IMD1);
				$rawFlag = 0;
				$zippedFlag = 0;
				foreach $file (@theFileNames){
					if($file eq "raw_data.bin"){
						$rawFlag = 1;}
					if($file eq "raw_data.bin.tar.gz"){
						$zippedFlag = 1;}
				}
				if($rawFlag && $zippedFlag){
					chdir($drc);
					print "Working in directory $drc \n";
					`rm raw_data.bin`;
					chdir($mainDir);
				}
				closedir(IMD1);
			}
		}
	}
	
	closedir(IMD);
}


