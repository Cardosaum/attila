#! /usr/perl/bin

use warnings;
use strict;
#~ ------------------------------------------------------------------------
#~ Programa: frequency_counter.pl
#~ Data: 24/04/2015
#~ Este programa recebe como entrada uma lista de sequ�ncias com 
#~ no seguinte formato: 
#~ >seq1
#~ seq
#~ #idcdr
#~ cdr*
#~ <
#~ Em que (<) marca o fim do arquivo. Aten��o, este programa l� as 
#~ strings considerando que estas est�o em apenas uma linha, n�o est�o 
#~ subdivididas em linhas de n caracteres. O programa armazena os dados
#~ em um vetor hashes. A chave do hash � a substring contendo cdrs, 
#~ e o campo freq armazena o n�mero de sequ�ncias �nicas igual a uma dada
#~ chave. O programa usa um outro hash para armazenar todas as sequ�ncias e
#~ seus respectivos ids que possuem a mesma subtring de cdrs. 
#~ Finalmente, � impressa uma lista ordenada por frequ�ncia, contendo 
#~ o id da cdr e sua frequ�ncia, a substring contendo cdrs, e todas as 
#~ sequ�ncias e seus ids que contem a substring.
#~ ------------------------------------------------------------------------

	#~ Passo. Declara��o de vari�veis
	my ($id,$seq,$idcdr,$cdr,$line,$temp);
	my ($flag,$len,$count,$c,$i,$k,$libsize,$freq);
	my ($input,$inputfiltered,$output) = @ARGV;
	my (@seqs);
	my %lista;	
	$cdr = "";
	$i = 0;
	$count = 0;
	$k = 0;
#~ -----------------------------------------------------------------------
								#~ Leitura
#~ -----------------------------------------------------------------------

	#~ Passo. Se o usu�rio n�o informou todos os argumentos
	if (@ARGV < 3){
		#~ Passo. Imprima mensagem de erro e pare a execu��o
		die "frequency_counter3.pl: perl frequency_counter3.pl <inputfile_protein.fasta> <input_filtered.fasta> <outputfile>\n";
	}
	#~ Passo. Abra o arquivo de entrada
	open IN, "<$input" or die "Could not find input file $input\n";

	#~ Passo. Enquanto houver entrada
	while($line = <IN>){
		#~ Passo. Remova o '\n' da linha
		chomp($line);
		#~ Passo. Se linha for do id
		if($line =~ /^>(.+)/ || $line =~ /</){
			#~ Passo. Se tiver terminado de ler os dados da sequ�ncia anterior
			if ($cdr =~ /(.+)\*$/){
				#~ Passo. Obtenha o tamanho da substring cdr sem "*"
				$len = length($cdr) - 1;
				#~ Passo. Copie $cdr sem "*" para $temp
				$temp = substr($cdr,0,$len);
				#~ Passo. Busque a cdr nos hashes do elemento atual
				$c = busca_cdr(\$temp,\$k,\@seqs);
				#~ Passo. Se a cdr foi encontrada
				if ($c != -1){
					#~ Passo. Incremente a frequ�ncia da cdr
					$seqs[$c]{$temp}{freq} = $seqs[$c]{$temp}{freq} + 1;
					#~ Passo. Armazene id e sequ�ncia que contem a substring de cdrs
					$lista{$temp}{id}[$seqs[$c]{$temp}{freq}-1] = $id;
					$lista{$temp}{seq}[$seqs[$c]{$temp}{freq}-1] = $seq;
				}
				else{
					#~ Passo. Inicialize a frequ�ncia da cdr
					$seqs[$i]{$temp}{freq} = 1;;
					#~ Passo. Armazene id e sequ�ncia que contem a substring de cdrs
					$lista{$temp}{id}[0] = $id;
					$lista{$temp}{seq}[0] = $seq;
					#~ Passo. Incremente o contador de sequ�ncias
					$count++;
				}				
				#~ Passo. Armazene o id da cdr em lista
				$lista{$temp}{idcdr} = $idcdr;
				#~ Passo. Se 10000 j� foram armazenadas no hash do elemento atual
				if ($count % 10000 == 0){
					#~ Passo. V� para o pr�ximo elemento do vetor de hashes
					$i++;
					#~ Passo. Volte um elemento
					$k = $i - 1;
				}
				else{
					#~ Passo. Copie $i para $k
					$k = $i;
				}				
			}
			#~ Passo. Armazene linha na tempor�ria id
			$id = $line;
			#~ Passo. Inicialize seq
			$seq = "";
			#~ Passo. Inicialize flag
			$flag = 1;
		}
		#~ Passo. Sen�o se a linha for do idcdr
		elsif($line =~ /^#(.+)/){
			#~ Passo. Armazene linha na tempor�ria idcdr
			$idcdr = $1;
			#~ Passo. Inicialize cdr
			$cdr = "";
			#~ Passo. Inicialize flag
			$flag = 2;
		}
		#~ Passo. Sen�o se a linha for de seq
		elsif($flag == 1){
			#~ Passo. Armazene linha em seq
			$seq = $line;
		}
		else{
			#~ Passo. Armazene a linha em cdr
			$cdr = $line;
		}
	}
	
	#~ Passo. Feche o arquivo de entrada
	close(IN);
#~ ------------------------------------------------------------------------
					#~ Processamento
#~ ------------------------------------------------------------------------	
	
	#~ Passo. Obtenha o tamanho da biblioteca de sequ�ncias filtradas
	$libsize = `grep -cP "^>" $inputfiltered`;
	chomp($libsize);
	$libsize = $libsize + 0;
	print "libsize : $libsize\n";
	#~ Passo. Se o �ltimo elemento do vetor totaliza 10000 hashes
	if ($count % 10000 == 0){
		#~ Passo. Volte um elemento
		$i = $i -1;
	}
	
	#~ Passo. Inicialize $c
	$c = 0;
	#~ Enquanto existirem elementos no vetor de hash
	while($c <= $i){
		#~ Passo. Para cada chave cdr
		foreach $line (keys %{$seqs[$c]}){
			#~ Passo. Copie a frequ�ncia da cdr para o hash lista
			$lista{$line}{freq} = $seqs[$c]{$line}{freq};
		}
		#~ Passo. V� para o pr�ximo elemento do vetor
		$c++;
	}
#~ ------------------------------------------------------------------------
						#~ Escrita
#~ ------------------------------------------------------------------------	
	#~ Passo. Abra o arquivo de sa�da
	open OUT, ">$output" or die "frequency_counter2.pl: Could not create output file $output\n";
	#~ Passo. Para cada cdr
	foreach $k (sort{$lista{$b}{freq} <=> $lista{$a}{freq}} keys %lista){
		#~ Passo. Inicialize o �ndice dos vetores de ids e seqs
		$c = 0;
		#~ Passo. Calcule a frequ�ncia relativa da cdr
		$freq = ($lista{$k}{freq} / $libsize) * 100000000000000;
		#~ Passo. Imprima o id da cdr e frequencia
		print OUT "#$lista{$k}{idcdr}|$libsize|$freq\n";
		#~ Passo. Imprima cdr 
		print OUT "*$k\n";
		#~ Passo. Enquanto existirem ids e sequ�ncias associadas a cdr atual
		while($c < $lista{$k}{freq}){
			#~ Passo. Imprima o id da sequ�ncia
			print OUT "$lista{$k}{id}[$c]\n";
			#~ Passo. Imprima a sequ�ncia
			print OUT "$lista{$k}{seq}[$c]\n";
			#~ Passo. V� para o proximo id e pr�xima sequ�ncia
			$c++;
		}
	}
	#~ Passo. Feche o arquivo de sa�da
	close(OUT);
	
exit(0);

#~ ------------------------------------------------------------------------
				#~ Defini��o de subrotinas
#~ ------------------------------------------------------------------------

sub busca_cdr{
	
	#~ Passo. Declara��o e atribui��o de vari�veis

	my $cdr = shift or die "frequency_counter2.pl: Could not receive string reference\n";
	$cdr = ${$cdr};
	my $i = shift or die "frequency_counter2.pl: Could not receive integer index reference\n";
	$i = ${$i};
	my $temp = shift or die "frequency_counter2.pl: Could not receive array\n";
	my @seqs = @{$temp};
	my $flag = -1;
	my $c = 0;
	#~ Passo. Enquanto n�o chegar ao fim do vetor e n�o achar a cdr
	while($c <= $i && $flag == -1)
	{
		#~ Passo. Se a cdr for chave de algum hash do elemento atual
		if(exists $seqs[$c]{$cdr})
		{
			#~ Passo. Modifique flag
			$flag = $c;
		}
		else
		{
			#~ Passo. V� para o pr�ximo elemento do vetor de hashes
			$c++;
		}
	}
	
	#~ Passo. Retorne o �ndice $c
	return $flag;
}
		
