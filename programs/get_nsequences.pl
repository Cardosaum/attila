#! usr/bin/perl

use warnings;
use strict;

#~ ----------------------------------------------------------------------
# Programa: get_nsequences.pl
#~ Este programa copia as primeiras n sequ�ncias do arquivo de entrada, e
#~ as escreve num arquivo de sa�da.
#~ ----------------------------------------------------------------------

	#~ Passo. Declara��o de vari�veis
	my ($number,$line,$header,%seq,$n,$seq);
	my ($input,$output,$total);
	$number = 0;
	$header = "";
	$seq = "";

	# Passo. Leitura

	if (@ARGV < 3){	
		die "get_nsequences: perl get_nsequences.pl <inputfile.fasta> <outputfile> <integer>\n";
	}

	($input,$output,$n) = @ARGV;
		
	# Passo. Abra o arquivo de entrada
	open IN, "<$input" or die "get_nsequences: Could not open input file $input: $! \n";

	#~ Passo. Conte o n�mero de sequ�ncias do arquivo
	$total = `grep -c -P ^'>' $input`;
	
	#~ Passo. Se o n�mero de sequ�ncias solicitado for maior que o total de sequ�ncias
	if($n > $total){
		#~ Passo. Atribua o conte�do de $total a $n
		$n = $total;
	}
	#~ Passo. Enquanto n�o ler n sequ�ncias
	while ($n > 0){
		#~ Passo. Leia a linha
		$line = <IN>;
		#~ Passo. Se tiver terminado de ler os dados da sequ�ncia anterior
		if(length($seq) > 1){
			# Passo. Se linha for do id ou tiver terminado de ler o arquivo
			if((!defined($line)) || ($line =~ /^>(.+)/)){
				#~ Passo. Armazene os dados da sequ�ncia anterior no hash
				$seq{$header}{read} = $seq;
				$seq{$header}{number} = $number;
				$number++;
				$n--;
			}
		}
		#~ Passo. Se a linha tiver algum conte�do definido
		if(defined($line)){
			#~ Passo. Remova o '\n' da linha
			chomp($line);
			#~ Passo. Se a linha for do id
			if($line =~ /^>(.+)/){
			#~ Passo. Armazene a linha em $header
			$header = $1;
			#~ Passo. Inicialize $seq
			$seq = "";			
			}
			else{
				# Passo. Concatene a linha em $seq
				$seq .= $line;
			}
		}
	}	
	close IN;

	open OUT, ">$output" or die "get_nsequences: Could not open output file $output: $!\n";

	
	foreach $line (sort{$seq{$a}{number} <=> $seq{$b}{number}} keys(%seq)){
		#~ Passo. Imprima o id e a sequ�ncia no arquivo de sa�da
		print OUT ">$line\n$seq{$line}{read}\n";
	}
	
	close OUT;
	
	exit(0);
