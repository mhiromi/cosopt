#!/usr/bin/perl
#
# 引数: １個目はソース、２個目はリザルト、３個目はイテレーション回数
#
print @ARGV;
my $source = (@ARGV) ? shift @ARGV : die "Input 1 argument at least";
my $result = (@ARGV) ? shift @ARGV : "result";
my $iteration = (@ARGV) ? shift @ARGV : 1;
my $sigma = (@ARGV) ? shift @ARGV : 0.3;
print $iteration . "\n";
$iteration = int($iteration) >= 1 ? int($iteration) : 1;
open(MPD_HOSTS,"$ENV{HOME}/mpd.hosts");
my $queues;
while(<MPD_HOSTS>){
	chomp;
	#$queues = $queues . $_ . ".q@" . $_;
	$queues = $queues . $_ . ".q";
	$queues = $queues . ",";
}
$queues =~ s/^(.*),$/$1/;

print $iteration . "\n";
for($i=1;$i <= $iteration; $i+=1){
$output = $result . "_" . $i . ".txt";
system("qsub -o $output -q $queues mpiscript.sh " . $source . " " . $sigma);
}
