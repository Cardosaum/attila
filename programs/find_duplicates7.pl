#! usr/bin/perl

use warnings;
use strict;

#----------------------------------------------------------------------------------------------------------------------
# Programa: find_duplicates6.pl
# Este programa recebe dois arquivos de entrada, ambos contendo uma lista
# de sequ�ncias de cdr com suas respectivas frequ�ncias e ids. O programa busca
# sequ�ncias cuja frequ�ncia aumentou do arquivo 1 em rela��o ao arquivo 2, e as
# imprime num arquivo de sa�da. Este programa volta a fazer uma coisa, que eu pensei que
# era gambiarra, mas em perl n�o �. Ent�o a estrutura de dados onde s�o armazenadas as
# sequ�ncias � um vetor de hash, por�m, as chaves dos hashes, diferente das 
# outras vers�es, s�o as pr�prias sequ�ncias das cdrs.Este programa imprime uma lista ordenada de
# sequ�ncias de cdrs de acordo com o fold change (freq2/freq1).
##----------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------------------
# 					In�cio da fun��o principal
#------------------------------------------------------------------------------------------------------------------------

	# Passo. Se o programa receber menos argumentos que o necess�rio
	if (@ARGV < 3)
	{
	# 	Passo. Pare a execu��o do programa e imprima mensagem de erro
		die "find_duplicates7.pl: perl find_duplicates7.pl <inputfile1> <inputfile2> <outputfile>\n";
	}

	# Passo. Declara��o de vari�veis 
	my ($input1,$input2,$output) = @ARGV;
	my ($line,$header,$freq,$number,$x,$y);
	my ($i,$c,$f,$counter,$libsize,$temp);
	my ($idseq,$seq,$cdr,$k,$frame);
	my @cdrs;
	my %allcdrs;

	# Passo. Inicialize o contador de sequ�ncias lidas
	$number = 0;
	# Passo. Inicialize o �ndice do vetor de hashes
	$i = 0;

#---------------------------------------------------------------------------------------------------------------------------
#					Leitura do arquivo 1
#---------------------------------------------------------------------------------------------------------------------------

	# Passo. Abra o arquivo de entrada
	open IN, "<$input1" or die "find_duplicates7.pl: Could not open input file $input1: $!\n";

	# Passo. Enquanto houver entrada
	while ($line = <IN>){
	# 	Passo. Remova o caracter \n do fim da linha
		chomp($line);

	# 	Passo.Se a linha come�ar com '#'
		if ($line =~ /^#(.+)\|FRAME:(\d)\|(.+)\|(\d+)\|((\d+)\.?(\d+))$/){		
	# 		Passo. Armazene a primeira express�o regular em $header
			$header = $1;
			#~ Passo. Armazene a segunda express�o regular como tamanho da biblioteca
			$libsize = $4;
	# 		Passo. Armazene a terceira express�o regular em $freq
			$freq = $5;
		}
		else{
	# 		Passo. Caso a linha comece com '*'
			if($line =~ /^\*(.+)/){
	# 			Passo. Armazene $freq no campo freq, usando a sequ�ncia como chave 
				$cdrs[$i]{$1}{freq1} = $freq;
	# 			Passo. Inicialize o campo freq2
				$cdrs[$i]{$1}{freq2} = 0;
	# 			Passo. Incremente o contador de Leitura
				$number++;
		# 		Passo. Se 10000 sequ�ncias j� foram lidas
				if ($number % 10000 == 0){
	#		 		Passo. V� para o pr�ximo hash do vetor
					$i++;
				}
			}

		}

	}

	# Passo. Feche o arquivo 1
	close IN;

#---------------------------------------------------------------------------------------------------------------------------
#					Leitura do arquivo 2 e processamento
#---------------------------------------------------------------------------------------------------------------------------
	# Passo. Abra o arquivo 2
	open IN, "<$input2" or die "find_duplicates7.pl: Could not open input file $input2: $!\n";

	$k = -1;
	# Passo. Enquanto houver entrada
	while ($line = <IN>){
	# 	Passo. Remova o caracter \n do fim da linha
		chomp($line);
		
		# 	Passo.Se a linha come�ar com '#'
		if ($line =~ /^#(.+)\|FRAME:(\d)\|(.+)\|(\d+)\|((\d+)\.?(\d+))$/){		
	# 		Passo. Armazene a quinta express�o regular em $freq
			$freq = $5;
		}
		#~ Passo. Sen�o se a linha for da susbtring de cdrs
		elsif ($line =~ /^\*(.+)/){
			#~ Passo. Armazene a linha em cdr
			$cdr = $1;
	# 		Passo. Se $number � m�ltiplo de 10000
			if ($number % 10000 == 0){
	# 			Passo. Atribua como �ltimo �ndice $i-1
				$f = $i - 1;
			}
			else{
	# 			Passo. Atribua como �ltimo �ndice $i
				$f = $i;
			}
	#		Passo. Busque a cdr2 no hash cdrs
			$c = busca_cdr(\@cdrs,$1,\$f) ;
	# 					
	# 		Passo. Se a sequ�ncia do ciclo 3 existir no ciclo 1
			if($c != -1){
				#~ Passo. k recebe o �ndice da cdr
 				$k = $c;
			}
			else{
				#~ Passo. k recebe o �ndice atual do vetor de hashes
				$k = $i;
				#~ Passo. Armazene 1/libsize no campo freq1 usando a sequ�ncia como chave
				$cdrs[$k]{$1}{freq1} = (1 / $libsize) * 100000000000000 ;
#				Passo. Incremente o contador de leitura
				$number++;
			}
			
			#~ Passo. Armazene freq no campo freq2
			$cdrs[$k]{$1}{freq2} = $freq;
# 			Passo. Se a freq2 for maior que a frequ�ncia 1
			if ($cdrs[$k]{$1}{freq2} > $cdrs[$k]{$1}{freq1}){
# 				Passo. Armazene o fold change da sequ�ncia no hash dif
				$temp = $cdrs[$k]{$1}{freq2} / $cdrs[$k]{$1}{freq1};
				$allcdrs{$1}{dif} = sprintf("%.4f", $temp);
# 				Passo. Armazene o �ndice do vetor de hashes onde se encontra esta sequ�ncia
				$allcdrs{$1}{index} = $k;
				#~ Passo. Inicialize os campo id e seq de allcdrs
				$allcdrs{$1}{id} = "";
				$allcdrs{$1}{seq} = "";
				
			}
# 			Passo. Se $number � m�ltiplo de 10000
			if ($number % 10000 == 0){
#	 				Passo. V� para o proximo elemento do vetor
					$i++;
			}
		}
		#~ Passo. Se a linha for do id de uma sequ�ncia completa
		elsif($line =~ /^>(.+)/){
				#~ Passo. Armazene linha em idseq
				$idseq = $1;
		}
		else{
			#~ Passo. Armazene linha em seq
			$seq = $line;
			if(exists $allcdrs{$cdr}){
				if(length($seq) > length($allcdrs{$cdr}{seq})){
					$allcdrs{$cdr}{id} = $idseq;
					$allcdrs{$cdr}{seq} = $seq;
				}
			}
		}
	}
	# Passo. Feche o arquivo 2
	close(IN);

	# Passo. Se $number � m�ltiplo de 10000
	if ($number % 10000 == 0){
	# 	Passo. Atribua $i -1 como �ltimo �ndice do vetor
		$f = $i - 1;
	}
	else{
	# 	Passo. Atribua $i como �ltimo �ndice do vetor
		$f = $i;
	}

#---------------------------------------------------------------------------------------------------------------------------
#						Escrita
#---------------------------------------------------------------------------------------------------------------------------

	# Passo. Abra o arquivo de sa�da
	open OUT, ">$output" or die "find_duplicates7.pl: Could not open output file $output\n";

	# Passo. Obtenha a lista ordenada de chaves do hash %allcdrs de acordo com o valor do hash dif 
	foreach $header (sort{$allcdrs{$b}{dif} <=> $allcdrs{$a}{dif}} keys %allcdrs){
	# 	Passo. Atribua o �ndice do vetor de hash � variavel tempor�ria $c
		$c = $allcdrs{$header}{index};
		#~ $cdrs[$c]{$header}{freq1} = $cdrs[$c]{$header}{freq1} / 100000000000000 ;
		#~ $cdrs[$c]{$header}{freq2} = $cdrs[$c]{$header}{freq2} / 100000000000000 ;
	# 	Passo. Imprima o id, a frequ�ncia do ciclo 1, a frequ�ncia do ciclo 3 e o fold change da sequ�ncia atual
		print OUT ">$allcdrs{$header}{id}|FOLD-CHANGE:$allcdrs{$header}{dif}|P1:$cdrs[$c]{$header}{freq1}|P2:$cdrs[$c]{$header}{freq2}\n";
	# 	Passo. Imprima a sequ�ncia da cdrs
		print OUT "$allcdrs{$header}{seq}\n";
	}

	# Passo. Feche o arquivo de sa�da
	close OUT;

	exit (0);

#----------------------------------------------------------------------------------------------------------------------------------------
# Subrotina busca_cdr: Esta subrotina recebe um vetor de hashes, uma string, e o �ltimo �ndice cujo elemento foi preenchido
# no vetor. A subrotina busca no vetor de hashes uma chave que seja igual a string recebida. Se achar uma chave, uma vari�vel
# auxiliar armazena o valor do �ndice do vetor onde est� o hash com tal chave. A subrotina devolve a vari�vel auxiliar, que
# ter� valor igual a -1 caso n�o seja encontrada nenhuma chave igual � string, ou o valor do �ndice do vetor onde foi encontrada a chave.
#----------------------------------------------------------------------------------------------------------------------------------------
sub busca_cdr{

	# Passo. Receba as vari�veis da main
	my $temp = shift or die "find_duplicates7.pl: Could not receive array of hashes\n";
	my @cdrs = @$temp;
	my $cdr2 = shift or die "find_duplicates7.pl: Could not receive string\n";
	$temp = shift  or die "find_duplicates7.pl: Could not receive last array index";
	my $f = $$temp;
	my ($c,$flag);

	# Passo. Inicialize flag
	$flag = -1;
	# Passo. Inicialize $c
	$c = 0;
	# Passo. Enquanto n�o achar a chave em um hash e enquanto n�o chegar ao fim do vetor
	while($flag == -1 && $c <= $f){
	# 	Passo. Se a chave existe no hash atual
		if (exists $cdrs[$c]{$cdr2}){
	# 		Passo. Altere flag
			$flag = $c;
		}
		else{
	# 		Passo. V� para o pr�ximo elemento do vetor
			$c++;
		}
	}

	# Passo. Retorne o �ndice encontrado
	return $flag;

}


