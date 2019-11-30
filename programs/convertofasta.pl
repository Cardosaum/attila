#! usr/bin/perl

use warnings;
use strict;
#--------------------------------------------------------------------------------------------------------------
# Programa:convertofasta.pl
# Este programa recebe como entrada um arquivo que contem uma lista de sequ�ncias dotadas
# de id e sequ�ncia. A sequ�ncia se encontra em formato de colunas numeradas, em que cada
# res�duo de amino�cido recebe um n�mero de acordo com o esquema de numera��o de imunoglobulinas
# de Kabat. O programa concatena cada um dos caracteres que representam os amino�cidos para formar
# uma sequ�ncia de dom�nio vari�vel completa e imprime a sequ�ncia em formato fasta.
#--------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------
#				In�cio da fun��o principal
#--------------------------------------------------------------------------------------------------------------

# Passo. Declara��o de vari�veis
my ($line,$number,$aux1);
my %cdrs;
# Passo. Se o usu�rio entrar com menos argumentos que o necess�rio
if (@ARGV < 2)
{
# 	Passo. Imprima mensagem de erro e pare a execu��o do programa
	die "convertfasta.pl: perl convertofasta.pl <inputfile> <outputfile> \n";
}
my ($input,$output) = @ARGV;

#--------------------------------------------------------------------------------------------------------------
#					Leitura
#--------------------------------------------------------------------------------------------------------------

# Passo. Abra o arquivo de entrada
open IN, "<$input" or die "convertfasta.pl: Could not open input file $input:$!\n";

# Passo. Inicialize o contador da ordem de leitura
$number = 0;

# Passo. Enquanto houver entrada
while ($line = <IN>)
{
# 	Passo. Remova o "\n" da linha
	chomp($line);
# 	Passo. Se a linha for do id da sequ�ncia
	if ($line =~ /^>(.+)/)
	{
# 		Passo. Armazene o id em $aux1
		$aux1 = $1;
# 		Passo. Inicialize o hash sequ�ncia usando o id como chave
		$cdrs{$aux1}{seq} = "";
# 		Passo. Atribua o contador de leitura para o hash number
		$cdrs{$aux1}{number} = $number;
# 		Passo. Incremente o contador de leitura
		$number++;
	}
# 	Passo. Sen�o
	else
	{
# 		Passo. Se a linha for de uma coluna de res�duo numerado
		if ($line =~ /([A-Z]?)$/)
		{
# 			Passo. Concatene o �ltimo caracter da linha no campo seq
			if(defined($1)){
				$cdrs{$aux1}{seq} .= $1;
			}
		}
	}
}

# Passo. Feche o arquivo de sa�da
close (IN);

#--------------------------------------------------------------------------------------------------------------
#					Escrita
#--------------------------------------------------------------------------------------------------------------
# Passo. Abra o arquivo de sa�da
open OUT, ">$output" or die "convertfasta.pl: Could not open output file $output:$!\n";

# Para cada id do hash
foreach $line (sort{$cdrs{$a}{number} <=> $cdrs{$b}{number}} keys %cdrs)
{
# 	Passo. Imprima o id
	print OUT ">$line\n";
# 	Passo. Imprima a sequ�ncia
	print OUT "$cdrs{$line}{seq}\n";
}

# Passo. Feche o arquivo
close (OUT);

exit(0);

