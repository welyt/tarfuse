#!/usr/bin/env perl -w
use strict;
use Getopt::Long;
use Cwd;
use File::Basename;
use File::Spec;
use configdata;

my @usage;
push @usage, "Usage: ".basename($0)." [options]\n";
push @usage, "blast the result of fusion.out of tophatfusion.\n";
push @usage, "  -h, --help      Displays this information\n";
push @usage, "  -c, --config    Configuration Filename\n";
push @usage, "  -t, --tophat    the output dir of tophatfusion\n";
push @usage, "  -o, --outdir    the output directory name\n";
push @usage, "  -1, --fq1       the first fastq of pair end\n ";
push @usage, "  -2, --fq2       the second fastq of pair end\n";
push @usage, "  -p, --prefix    the output file prefix name\n";


my $help;
my $config_file;
my $tophat;
my $outdir;
my $fq1;
my $fq2;
my $prefix;

GetOptions
(
	'help'	      => \$help,
	'config=s'    =>\$config_file,
	'tophat=s'    => \$tophat,
	'outdir=s'    => \$outdir,
	'1=s'       => \$fq1,
	'2=s'       => \$fq2,
	'prefix=s'    => \$prefix,
);

not defined $help  or die @usage;
defined $tophat    or die @usage;
defined $outdir    or die @usage;
defined $fq1       or die @usage;
defined $fq2       or die @usage;

mkdir $outdir if not -d $outdir;
-e $config_file or die "Error: Unable to find config file $config_file\n";
if (not defined $prefix){
	$prefix = "";
	my @output_splitdir = File::Spec->splitdir($outdir);
	while ($prefix eq "" and scalar @output_splitdir > 0){
		$prefix = pop(@output_splitdir);
		defined $prefix or die "Error: Unable to infer library name from output director $outdir\n";
	}
}


my $config=configdata->new();
$config->read($config_file);
my $genome_version          = $config->get_value("genome_version");
my $genome_fasta            = $config->get_value("genome_fasta");
my $nt_database             = $config->get_value("nt_database");
my $human_genomic_database  = $config->get_value("human_genomic_database");
my $reads_span_number       = $config->get_value("reads_span_number");
my $left_reads_cov_length   = $config->get_value("left_reads_cov_length");
my $right_reads_cov_length  = $config->get_value("right_reads_cov_length");
my $nt_score                = $config->get_value("blast_nt_score");
my $hg_score                = $config->get_value("blast_HG_score");
my $two_score               = $config->get_value("blast_two_score");
my $num_threads             = $config->get_value("num_threads");


my $log           = "$outdir/$prefix\_log.out";
my $tophat_fusion = "$tophat/fusions.out";
open(FU,">$outdir/$prefix\_fusion_filter.out") or die $!;
my ($fusion100bp,$fusionRef,$fusioninfo) = getinfo($tophat_fusion,$reads_span_number,$left_reads_cov_length,$right_reads_cov_length,$log);
my ($blast2score)                        = blast($fusion100bp,$nt_database,$human_genomic_database,$num_threads,$prefix,$outdir,$log);

print FU "chr\tleft_chr_fusion_point\tright_chr_fusion_point\torientation\tspan_reads\tpair\tspan_pair\tcontradict_reads\tleft_cov_length\tright_cov_length\tdistribute\tleftseq_fusion\trightseq_fusion\tleft_cov\tright_cov\tblast_two_sseqid\tblast_two_length_pident_sum\tblast_nt_sseqid\tblast_nt_length_pident_sum\tblast_HG_sseqid\tblast_HG_length_pident_sum\n";

foreach my $id (sort {$a cmp $b} keys %{$fusioninfo}){
	if(exists $blast2score->{$id}){
		print FU $fusioninfo->{$id},"\t",$blast2score->{$id}->{'blastTwo'},"\t",$blast2score->{$id}->{'blastNt'},"\t",$blast2score->{$id}->{'blastHG'},"\n";
	}
}
my $tophat_bam = "$tophat/accepted_hits.bam";

fusionBam($tophat_bam,$fq1,$fq2,$fusionRef,$genome_fasta,$num_threads,$prefix,$outdir,$log);


sub getinfo{
	my($tophat_fusion,$reads_span_number,$left_reads_cov_length,$right_reads_cov_length,$log)=@_;

	open(LOG,">>$log") or  die $!;
	print LOG "EXTRACT information and simple filter from the tophatfusion out\n";
	my $fusion100bp;
	my $fusionRef;
	my $fusioninfo;
	
	open(FUSION,"$tophat_fusion") or die "Please input your fusion.out file of tophatfusion\n";
	while(my $line=<FUSION>){
		chomp $line;
		$line=~s/\@\t//g;
		my @cols=split(/\t+/,$line);
		next if($#cols<15);
		if($cols[4] >= $reads_span_number and $cols[8] >= $left_reads_cov_length and $cols[9] >= $right_reads_cov_length){
			my $id="$cols[0]\_$cols[1]\_$cols[2]";
			@{$fusion100bp->{$id}->{'left'}}  = split(/\s+/,$cols[12]);
			@{$fusion100bp->{$id}->{'right'}} = split(/\s+/,$cols[13]);
			my @chr=split(/\-/,$cols[0]);
			$fusionRef->{$id}->{'leftchr'}     = $chr[0];
			$fusionRef->{$id}->{'rightchr'}    = $chr[1];
			$fusionRef->{$id}->{'leftPoint'}   = $cols[1];
			$fusionRef->{$id}->{'rightPoint'}  = $cols[2];
			$fusionRef->{$id}->{'orientation'} = $cols[3];
			$fusioninfo->{$id}=join("\t",@cols[0 .. 10],@cols[12 .. 15]);
			
		}
	}
	print LOG "FINISH extrate information and simple filter\n";
	close(LOG);
	return($fusion100bp,$fusionRef,$fusioninfo);
}


sub blast{
	
	my ($fusion100bp,$nt_database,$HG_database,$num_threads,$prefix,$outdir,$log)=@_;
	open(LOG,">>$log") or die $!;
	print LOG "BEGIN blast nt and human genomic\n";
	my $blast2score;
	foreach my $id (sort {$a cmp $b} keys %{$fusion100bp}){
		my $leftseq=join('',@{$fusion100bp->{$id}->{'left'}});
		my $rightseq=join('',@{$fusion100bp->{$id}->{'right'}});
		my $fusionseq=join('',$fusion100bp->{$id}->{'left'}->[0],$fusion100bp->{$id}->{'right'}->[1]);
		#my $dir=dirname(Cwd::abs_path($outdir));
		my $tmp="$outdir/$prefix\_tmp";
		mkdir("$tmp");
		my $left_seq_file    = "$tmp/$id\_leftseq";
		my $right_seq_file   = "$tmp/$id\_rightseq";
		my $fusion_seq_file  = "$tmp/$id\_fusionseq";
		open(LEFTTEMP,">$left_seq_file") or die $!;
		print LEFTTEMP ">$id\_left\n$leftseq\n";
		open(RIGHTTEMP,">$right_seq_file") or die $!;
		print RIGHTTEMP ">$id\_right\n$rightseq\n";
		open(FUSIONTEMP,">$fusion_seq_file") or die $!;
		print FUSIONTEMP ">$id\_100bp\_fusion\n$fusionseq\n";
		`makeblastdb -in $tmp/$id\_rightseq -dbtype nucl`;
		print LOG "BLASTING $id  ";
		my @blastTwo=`blastn -query $left_seq_file -subject $right_seq_file -outfmt \"6 qseqid sseqid length pident\" -word_size 28`;
		my @blastNt =`blastn -query $fusion_seq_file -db $nt_database -outfmt \"6 qseqid sseqid length pident\" -word_size 28 -num_threads $num_threads`;
		my @blastHG =`blastn -query $fusion_seq_file -db $HG_database -outfmt \"6 qseqid sseqid length pident\" -word_size 28 -num_threads $num_threads`;
		my $blast_two   = blastFilter(\@blastTwo,$two_score);
		my $blast_nt    = blastFilter(\@blastNt,$nt_score);
		my $blast_HG    = blastFilter(\@blastHG,$nt_score);
		if($blast_two eq 'undef' or $blast_nt eq 'undef' or $blast_HG eq 'undef'){
			next;
		}else{
			$blast2score->{$id}->{'blastTwo'} = $blast_two;
			$blast2score->{$id}->{'blastHG'}  = $blast_HG;
			$blast2score->{$id}->{'blastNt'}  = $blast_nt;
		}

		
	}
	print LOG "FINISH blast\n";
	close(LOG);
	return($blast2score);
}


sub blastFilter{
	my ($blast_result,$blast_score)=@_;
	my $blast_line;
	if(@{$blast_result}){
		my @blast_best=split(/\s+/,$blast_result->[0]);
		my $blast_score_sum=$blast_best[2]+$blast_best[3];
		if($blast_score_sum < $blast_score){
			$blast_line=join("\t",$blast_best[1],$blast_score_sum);
		}else{
			$blast_line='undef';
		}
	}else{
		my @noblast = ('-') x 2;
		$blast_line=join("\t",@noblast);

	}
	return($blast_line);
}


sub fusionBam{
	my ($alignBam,$fq1,$fq2,$fusionRef,$genome_fasta,$num_threads,$prefix,$outdir,$log)=@_;
	open(LOG,">>$log") or die $!;
	print LOG "GENERATE fusion Bam file\n";
	my $fusionfq_id;	
	open(ALIG,"samtools view $alignBam |") or die $!;
	while(my $line=<ALIG>){
		chomp $line;
		if($line =~/XF\:/){
			my @cols=split(/\t+/,$line);
			$fusionfq_id->{$cols[0]}=1;
		}
	}
	my $fq1handle;
	my $fq2handle;
	if($fq1 =~/\.gz/ and $fq2 =~/\.gz/){
		$fq1handle="zcat $fq1 \|";
		$fq2handle="zcat $fq2 \|";
	}else{
		$fq1handle=$fq1;
		$fq2handle=$fq2;
	}
	open(IN1,"$fq1handle") or die $!;
	open(IN2,"$fq2handle") or die $!;
	open(FQ1,">$outdir/$prefix\_fusion_R1.fq") or die $!;
	open(FQ2,">$outdir/$prefix\_fusion_R2.fq") or die $!;

	while(1){
		my $readid1 = <IN1>;
		my $sequence1 = <IN1>;
		my $comment1 = <IN1>;
		my $quality1 = <IN1>;
		last if not defined $quality1;
		chomp($readid1);
		chomp($sequence1);
		chomp($comment1);
		chomp($quality1);
		my $readid2 = <IN2>;
		my $sequence2 = <IN2>;
		my $comment2 = <IN2>;
		my $quality2 = <IN2>;
		last if not defined $quality2;
		chomp($readid2);
		chomp($sequence2);
		chomp($comment2);
		chomp($quality2);
		my $readid=$readid1;
		$readid=~s/^@//;
		my @fqId=split(/\s+/,$readid);
		if(exists $fusionfq_id -> {$fqId[0]}){
			print FQ1 "$readid1\n$sequence1\n$comment1\n$quality1\n";
			print FQ2 "$readid2\n$sequence2\n$comment2\n$quality2\n";
		}
	}
	print LOG "GENERATE fusion genome for mapping\n";
	my($fusionG)=fusionGenome($fusionRef,$genome_fasta);
	open(FUG,">$outdir/$prefix\_fusion_genome.fa") or die $!;
	foreach my $id2(sort {$a cmp $b} keys %{$fusionG}){
		print FUG $fusionG->{$id2}->[0],$fusionG->{$id2}->[1];
	}
	print LOG "FINISH generate fusion genome\n";
	my $index       =`bwa index $outdir/$prefix\_fusion_genome.fa`;
	print LOG "FINISH index fusion genome\n $index\n";
	my $mapping     =`bwa mem -t $num_threads $outdir/$prefix\_fusion_genome.fa $outdir/$prefix\_fusion_R1.fq $outdir/$prefix\_fusion_R2.fq > $outdir/$prefix\_fusion_align.sam`;
	print LOG "FINISH mapping to fusion genome\n $mapping\n";
	my $sam2bam     =`samtools view -Sb $outdir/$prefix\_fusion_align.sam -o $outdir/$prefix\_fusion_align.bam`;
	my $sortbam     =`samtools sort $outdir/$prefix\_fusion_align.bam $outdir/$prefix\_fusion_index`;
	my $indexbam    =`samtools index $outdir/$prefix\_fusion_index.bam`;
	print LOG "FINISH index bam file \n $sam2bam\n$sortbam\n$indexbam\n";
	my $indexGenome =`samtools faidx $outdir/$prefix\_fusion_genome.fa`;
	print LOG "FINISH generate fusion bam file\n";
	close(LOG);
}


sub fusionGenome{
	my ($fusionRef,$genomefa)=@_;
#my $genomefa="$genomedir/WholeGenomeFasta/genome.fa";
	my $fusionG;
	foreach my $id (sort {$a cmp $b} keys %{$fusionRef}){	
		if($fusionRef->{$id}->{'orientation'} eq 'ff'){
			$fusionRef->{$id}->{'leftStart'}=$fusionRef->{$id}->{'leftPoint'}-998;
			$fusionRef->{$id}->{'leftEnd'}=$fusionRef->{$id}->{'leftPoint'}+1;
			open(LEFTFA,"samtools faidx $genomefa $fusionRef->{$id}->{'leftchr'}:$fusionRef->{$id}->{'leftStart'}-$fusionRef->{$id}->{'leftEnd'}|") or die $!;
			my @leftfa=<LEFTFA>;
			chomp(@leftfa);
			my $leftseq=join('',@leftfa[1 .. $#leftfa]);
			$fusionRef->{$id}->{'rightStart'}=$fusionRef->{$id}->{'rightPoint'}+1;
			$fusionRef->{$id}->{'rightEnd'}=$fusionRef->{$id}->{'rightPoint'}+1000;
			open(RIGHTFA,"samtools faidx $genomefa $fusionRef->{$id}->{'rightchr'}:$fusionRef->{$id}->{'rightStart'}-$fusionRef->{$id}->{'rightEnd'}|") or die $!;
			my @rightfa=<RIGHTFA>;
			chomp(@rightfa);
			my $rightseq=join('',@rightfa[1 .. $#rightfa]);
			$fusionG->{$id}->[0]=">$id\_ff\n";
			$fusionG->{$id}->[1]="$leftseq$rightseq\n";
		}
		if($fusionRef->{$id}->{'orientation'} eq 'fr'){
			$fusionRef->{$id}->{'leftStart'}=$fusionRef->{$id}->{'leftPoint'}-998;
			$fusionRef->{$id}->{'leftEnd'}=$fusionRef->{$id}->{'leftPoint'}+1;
			open(LEFTFA,"samtools faidx $genomefa $fusionRef->{$id}->{'leftchr'}:$fusionRef->{$id}->{'leftStart'}-$fusionRef->{$id}->{'leftEnd'}|") or die $!;
			my @leftfa=<LEFTFA>;
			chomp(@leftfa);
			my $leftseq=join('',@leftfa[1 .. $#leftfa]);
			$fusionRef->{$id}->{'rightStart'}=$fusionRef->{$id}->{'rightPoint'}-998;
			$fusionRef->{$id}->{'rightEnd'}=$fusionRef->{$id}->{'rightPoint'}+1;
			open(RIGHTFA,"samtools faidx $genomefa $fusionRef->{$id}->{'rightchr'}:$fusionRef->{$id}->{'rightStart'}-$fusionRef->{$id}->{'rightEnd'}|") or die $!;
			my @rightfa=<RIGHTFA>;
			chomp(@rightfa);
			my ($revComfa)=revCom(\@rightfa);
			$fusionG->{$id}->[0]=">$id\_fr\n";
			$fusionG->{$id}->[1]="$leftseq$revComfa->[1]\n";
		}
		if($fusionRef->{$id}->{'orientation'} eq 'rf'){
			$fusionRef->{$id}->{'leftStart'}=$fusionRef->{$id}->{'leftPoint'}+1;
			$fusionRef->{$id}->{'leftEnd'}=$fusionRef->{$id}->{'leftPoint'}+1000;
			open(LEFTFA,"samtools faidx $genomefa $fusionRef->{$id}->{'leftchr'}:$fusionRef->{$id}->{'leftStart'}-$fusionRef->{$id}->{'leftEnd'}|") or die $!;
			my @leftfa=<LEFTFA>;
			chomp(@leftfa);
			my ($revComfa)=revCom(\@leftfa);
			$fusionRef->{$id}->{'rightStart'}=$fusionRef->{$id}->{'rightPoint'}+1;
			$fusionRef->{$id}->{'rightEnd'}=$fusionRef->{$id}->{'rightPoint'}+1000;
			open(RIGHTFA,"samtools faidx $genomefa $fusionRef->{$id}->{'rightchr'}:$fusionRef->{$id}->{'rightStart'}-$fusionRef->{$id}->{'rightEnd'}|") or die $!;
			my @rightfa=<RIGHTFA>;
			chomp(@rightfa);
			my $rightseq=join('',@rightfa[1 .. $#rightfa]);
			$fusionG->{$id}->[0]=">$id\_rf\n";
			$fusionG->{$id}->[1]="$revComfa->[1]$rightseq\n";
		}
		if($fusionRef->{$id}->{'orientation'} eq 'rr'){
			$fusionRef->{$id}->{'leftStart'}=$fusionRef->{$id}->{'leftPoint'}+1;
			$fusionRef->{$id}->{'leftEnd'}=$fusionRef->{$id}->{'leftPoint'}+1000;
			open(LEFTFA,"samtools faidx $genomefa $fusionRef->{$id}->{'leftchr'}:$fusionRef->{$id}->{'leftStart'}-$fusionRef->{$id}->{'leftEnd'}|") or die $!;
			my @leftfa=<LEFTFA>;
			chomp(@leftfa);
			my ($revComfaLeft)=revCom(\@leftfa);
			$fusionRef->{$id}->{'rightStart'}=$fusionRef->{$id}->{'rightPoint'}-998;
			$fusionRef->{$id}->{'rightEnd'}=$fusionRef->{$id}->{'rightPoint'}+1;
			open(RIGHTFA,"samtools faidx $genomefa $fusionRef->{$id}->{'rightchr'}:$fusionRef->{$id}->{'rightStart'}-$fusionRef->{$id}->{'rightEnd'}|") or die $!;
			my @rightfa=<RIGHTFA>;
			chomp(@rightfa);
			my ($revComfaRight)=revCom(\@rightfa);
			$fusionG->{$id}->[0]=">$id\_rr\n";
			$fusionG->{$id}->[1]="$revComfaLeft->[1]$revComfaRight->[1]\n";
		}
	}
	return($fusionG);

}


sub revCom{
	my ($fa)=@_;
	my $revComfa;
	$revComfa->[0]=$fa->[0];
	my $seq=join('',@{$fa}[1 .. $#{$fa}]);
	$seq =~ tr/ATGCatgcNn/TACGtacgNn/;
	$seq=scalar reverse $seq;
	$revComfa->[1]=$seq;
	return($revComfa);
}
