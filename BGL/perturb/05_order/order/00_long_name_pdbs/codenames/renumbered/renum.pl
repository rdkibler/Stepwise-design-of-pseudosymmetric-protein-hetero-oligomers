#!/usr/bin/perl
sub SR{my $str = $_[0]; $str =~ s/ //g; return $str;}

my $n = 0;
my $old_count;
if($ARGV[0]){$n = $ARGV[0];}
# ATOM   9953 3HG2 THR A 202     -13.364  17.786  10.102  1.00  0.00           H
my $k = 1;
while($line = <STDIN>)
{
	chomp($line);
	my $record = SR(substr($line,0,6));
	if(substr($line,0,6) eq "HETATM" and substr($line,17,3) eq "MSE")
	{
		$line =~ s/^HETATM/ATOM  /;
		substr($line,17,3,"MET");
		substr($line,12,2,"SD") if substr($line,12,2) eq "SE";
	}
	if($record eq "ATOM")
	{
		my $chain = substr($line,21,1);
		my $atom = substr($line,12,4);
		my $resi = substr($line,17,3);
		my $resn = SR(substr($line,22,5));

		$new_count = $resn;
		if($old_count eq "" or $new_count != $old_count){$n++;$old_count = $new_count;}


		my $x = SR(substr($line,30,8));
		my $y = SR(substr($line,38,8));
		my $z = SR(substr($line,46,8));

		#ATOM      1  N   MET A   1     115.750  10.284  -2.481  1.00  0.00
		#ATOM    911 3HD2 LEU A  58     236.813-168.783  49.440  1.00  0.00
		printf "ATOM  %5d %s %3s %s%4d    %8.3f%8.3f%8.3f  %4.2f  %4.2f\n",$k,$atom,$resi,$chain,$n,$x,$y,$z,1,0;
		$k++;
	}
}
