#This program is a data miner written by Ioannis Kagiampakis. For bugs: ikagiampis@gmail.com Enjoy!
use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;
use GD::Simple;
use GD;

$| = 1;

#---Warning_log------------------------------------
local $SIG{__WARN__} = sub {
	my $message = shift;
	logger( 'warning', $message );
};

my $counter = 0;
count();

sub count {
	$counter = $counter + 42;
}

sub logger {
	my ( $level, $msg ) = @_;
	if ( open my $out, '>>', 'log.txt' ) {
		chomp $msg;
		print $out "$level - $msg\n";
	}
}

#--------------------------------------------------
sub main {
	my %opts;

	#d: directory
	#r: gene rank
	#t: T-test
	#p: P-value
	#q: Q-value
	#f: Fold Change
	#l: gene number

	getopts( 'd:r:t:p:q:f:l:', \%opts );

	if ( !checkusage( \%opts ) ) {
		usage();
		exit();
	}

	my $input_dir = $opts{"d"};

	my $gene_number = $opts{"l"};

	my @files = get_files($input_dir);

	my ( $full_list, $repeat_genes ) =
	  process_files( \@files, $input_dir, \%opts );

	#print Dumper ($repeat_genes);
	#print Dumper ($full_list);
	my $limit = find_limit_of_search( $repeat_genes, $gene_number );

	my $total_counter = compare_list( $full_list, $repeat_genes, $limit );

	create_wordcloud( $repeat_genes, $limit, $total_counter );
}
#----------------------------------------------
sub find_limit_of_search {
	my ( $repeat_genes, $gene_number ) = @_;

	my $limit = 0;

	if ( defined($gene_number) ) {

		my $counter = 1;

		foreach my $repeats (
			sort { $repeat_genes->{$b} <=> $repeat_genes->{$a} }
			keys %{$repeat_genes}
		  )
		{
			if ( $gene_number >= $counter ) {

				$limit = $repeat_genes->{$repeats};
			}
			$counter++;
		}
	}
	return $limit;
}
#----------------------------------------------
sub create_wordcloud {
	my ( $repeat_genes, $limit, $total_counter ) = @_;

	my $image_y = 800 + ( $total_counter * 5 );

	my $img = GD::Simple->new( 800, "$image_y" );

	my $counter_x     = 40;
	my $counter_y     = 60;
	my $cycle_counter = 0;
	my $color_count   = 0;
	my $color_loop    = 1;
	my @color         = ( "red", "blue", "green", "orange", "black" );
	my $color         = 'black';

	foreach my $repeats_keys (
		sort { $repeat_genes->{$b} <=> $repeat_genes->{$a} }
		keys %{$repeat_genes}
	  )
	{

		if ( $repeat_genes->{$repeats_keys} >= $limit ) {

			my $count    = $repeat_genes->{$repeats_keys};
			my $fontsize = 10 * $count;

			if ( $color_count > $count ) {
				$cycle_counter++;
				$color_count = $count;
				$color       = $color[$cycle_counter];
				$img->fgcolor("$color");
			}

			if ( $color_count == $count ) {
				$img->fgcolor("$color");
			}

			if ( $cycle_counter == 0 ) {
				$img->fgcolor('red');
				$color_count = $count;
			}

			if ( $cycle_counter == 5 ) {
				$img->fgcolor('black');
			}
			my $font = $img->font( 'C:\Windows\Fonts\arialbd.ttf', $fontsize );

			$img->moveTo( "$counter_x", "$counter_y" );

			#$img->fontMetrics("$fontsize","$repeats_keys");
			$img->string("$repeats_keys");

		   #print "$repeats_keys, $count, $counter_x, $counter_y, $font_size\n";

			my ( $x, $y ) = $img->curPos;
			$counter_x = $x + 10;
			if ( $counter_x >= 430 ) {

				$counter_y = $y + 20;
				$counter_x = 60;
			}

		}
	}

	open my $out, '>', 'img.png' or die;
	binmode $out;
	print $out $img->png;

}
#----------------------------------------------
sub compare_list {

	my ( $full_list, $repeat_genes, $limit ) = @_;

	my $total_counter = 0;

	open( MYFILE, '>>final_data.csv' )
	  or die "Could not open file final_data.csv";

	print MYFILE
"\"Gene_Symbol\",\"Cancer Type\",\"Gene Rank\",\"Gene Name\",\"Reporter ID\",\"T-test\",\"P-value\",\"Q-value\",\"Fold Change\",\"File Name.\",\n";

	foreach my $repeats_n (
		sort { $repeat_genes->{$b} <=> $repeat_genes->{$a} }
		keys %{$repeat_genes}
	  )
	{

		#print "$repeats_n"."->"."$repeat_genes->{$repeats_n}\n";

		if ( $repeat_genes->{$repeats_n} >= $limit ) {

			foreach my $file_name ( keys %{$full_list} ) {

				my $file = $full_list->{$file_name};

				if ( exists $file->{$repeats_n} ) {
					my $data = $file->{$repeats_n};
					print MYFILE "\""
					  . $data->{gene_symbol} . "\",\""
					  . $data->{cancer} . "\",\""
					  . $data->{gene_rank} . "\",\""
					  . $data->{gene_name} . "\",\""
					  . $data->{reporter_id} . "\",\""
					  . $data->{ttest} . "\",\""
					  . $data->{pvalue} . "\",\""
					  . $data->{qvalue} . "\",\""
					  . $data->{fold_change} . "\",\""
					  . $file_name . "\",\n";

					$total_counter++;
				}

			}

		}

	}

	close(MYFILE);

	return $total_counter;
}

#----------------------------------------------
sub process_files {
	my ( $files, $input_dir, $opts ) = @_;

	my %repeat_list;
	my %all_list;
	my %main_counter;

	foreach my $file (@$files) {

		print "Processing $file in $input_dir ... \n";
		my %names;
		my $filepath = "$input_dir/$file";

		open( INPUTFILE, $filepath ) or die "Unable to open $filepath\n";
		my $counter = 0;
		my $cancer_name;
	  OPT1: while ( my $line = <INPUTFILE> ) {

			if ( $counter == 0 ) {

				chomp $line;
				if ( $line =~ /\"(.+)\"/ ) {
					$cancer_name = $1;
				}
			}
			$counter++;

			if ( $line =~
/^\"(\d+)\",\"(.+)\",\"(.+)\",\"(.+)\","(.+)\","(.+)\","(.+)\","(.+)\"$/
			  )
			{

				my $gene_rank_d   = $1;
				my $gene_name_d   = $2;
				my $gene_symbol_d = $3;
				my $reporter_id   = $4;
				my $ttest_d       = $5;
				my $pvalue_d      = $6;
				my $qvalue_d      = $7;
				my $fold_change   = $8;

				if ( !con_gene_rank( $gene_rank_d, $opts ) ) {
					next OPT1;
				}

				if ( !con_ttest( $ttest_d, $opts ) ) {
					next OPT1;
				}

				if ( !con_pvalue( $pvalue_d, $opts ) ) {
					next OPT1;
				}

				if ( !con_qvalue( $qvalue_d, $opts ) ) {
					next OPT1;
				}

				if ( !con_fold_change( $fold_change, $opts ) ) {
					next OPT1;
				}

				my %data;
				$data{"gene_rank"}   = $gene_rank_d;
				$data{"gene_name"}   = $gene_name_d;
				$data{"gene_symbol"} = $gene_symbol_d;
				$data{"reporter_id"} = $reporter_id;
				$data{"ttest"}       = $ttest_d;
				$data{"pvalue"}      = $pvalue_d;
				$data{"qvalue"}      = $qvalue_d;
				$data{"fold_change"} = $fold_change;
				$data{"cancer"}      = $cancer_name;
				$data{"file_name"}   = $file;

				$names{$gene_symbol_d} = \%data;

				#print"$gene_symbol_d\n";
				if ( exists $repeat_list{$gene_symbol_d} ) {

					$main_counter{$gene_symbol_d} =
					  $main_counter{$gene_symbol_d} + 1;

					$repeat_list{$gene_symbol_d} =
					  $main_counter{$gene_symbol_d};

				}
				else {
					$repeat_list{$gene_symbol_d}  = 1;
					$main_counter{$gene_symbol_d} = 1;

				}
			}

		}
		close(INPUTFILE);

		#$main_counter++;
		$all_list{$file} = \%names;
	}

	#print Dumper (\%repeat_list);
	#print Dumper (\%all_list);
	return \%all_list, \%repeat_list;
}

#----------------------------------------------
sub con_gene_rank {

	my ( $gene_rank_d, $opts ) = @_;

	if ( defined( $opts->{"r"} ) ) {

		if ( $gene_rank_d <= $opts->{"r"} ) {
			return 1;
		}
		else {
			return 0;
		}

	}
	else {
		return 1;
	}

}

#----------------------------------------------
sub con_ttest {

	my ( $ttest_d, $opts ) = @_;

	if ( defined( $opts->{"t"} ) ) {

		if ( $opts->{"t"} >= 0 ) {
			if ( $ttest_d >= $opts->{"t"} ) {
				return 1;
			}
		}

		elsif ( $opts->{"t"} < 0 ) {
			if ( $ttest_d <= $opts->{"t"} ) {
				return 1;
			}
		}

		else {
			return 0;
		}
	}
	else {
		return 1;
	}

}

#----------------------------------------------
sub con_pvalue {

	my ( $pvalue, $opts ) = @_;

	if ( defined( $opts->{"p"} ) ) {

		if ( $pvalue <= $opts->{"p"} ) {
			return 1;
		}
		else {
			return 0;
		}

	}
	else {
		return 1;
	}
}

#----------------------------------------------
sub con_qvalue {

	my ( $qvalue, $opts ) = @_;

	if ( defined( $opts->{"q"} ) ) {

		if ( $qvalue <= $opts->{"q"} ) {
			return 1;
		}
		else {
			return 0;
		}

	}
	else {
		return 1;
	}
}

#----------------------------------------------
sub con_fold_change {

	my ( $fold_change, $opts ) = @_;

	if ( defined( $opts->{"f"} ) ) {

		if ( $opts->{"f"} >= 0 ) {
			if ( $fold_change >= $opts->{"f"} ) {
				return 1;
			}
		}

		elsif ( $opts->{"f"} < 0 ) {
			if ( $fold_change <= $opts->{"f"} ) {
				return 1;
			}
		}
		else {
			return 0;
		}

	}

	else {
		return 1;
	}
}

#----------------------------------------------

sub get_files {

	my $input_dir = shift;

	unless ( opendir( INPUTDIR, $input_dir ) ) {
		die "\n Unable to open durectory '$input_dir'\n";
	}

	my @files = readdir(INPUTDIR);

	#print Dumper( \@files );

	closedir(INPUTDIR);

	@files = grep( /\.csv$/i, @files );

	return @files;
}

#----------------------------------------------
sub checkusage {
	my $opts = shift;

	#d: directory
	#r: gene rank
	#t: T-test
	#p: P-value
	#q: Q-value
	#f: Fold Change

	my $input_dir   = $opts->{"d"};
	my $gene_rank   = $opts->{"r"};
	my $t_test      = $opts->{"t"};
	my $p_values    = $opts->{"p"};
	my $q_value     = $opts->{"q"};
	my $fold_change = $opts->{"f"};

	unless ( defined($input_dir) ) {
		return 0;
	}

	return 1;
}

#----------------------------------------------
sub usage {

	print <<USAGE;
	
usage: perl miner.pl <options>
	-d <directory>	specify directory in which to find csv files.
	-r <Gene Rank>
	-t <T-test value> 
	-p <P values>
	-q <Q values>
	-f <Fold Change>
	-l <How long you want the genes list to be>
	

example usage:
	
	#If you want to find the 100 most common genes, that show 1X Fold Change, you write
	perl miner.pl -d folder_name -f 1 -l 100
	
USAGE

}

#----------------------------------------------

main();
