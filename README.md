#################################################################################################################################################
#																		#
#                              				sf723_descriptors.py									#
#																		#
#################################################################################################################################################

É necessário instalar algumas dependências. Com o Anaconda é mais fácil implementar:

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
| Dependência 		| Script Conda/Apt-get 				| Observação						|
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
| Pandas 		| conda install -c anaconda pandas 		|							|
---------------------------------------------------------------------------------------------------------------------------------
| Biopython 		| conda install -c anaconda biopython 		|							|
---------------------------------------------------------------------------------------------------------------------------------
| DSSP 			| sudo apt-get install -y dssp 			|							|
---------------------------------------------------------------------------------------------------------------------------------
| Java 			| sudo apt install default-jre 			|							|
---------------------------------------------------------------------------------------------------------------------------------
| Rdkit 		| conda install -c rdkit rdkit 			|							|
---------------------------------------------------------------------------------------------------------------------------------
| OpenBabel 2.4 	| conda install -c conda-forge openbabel	| É obrigatório que seja a versão 2.4 para criar os 	|
|			|						| atributos do NNScore 2.0, a versão mais recente do	|
|			|						| OpenBabel gera problemas porque mudou completamente  	|
|			|						| a maneira de ler as moléculas.			|
---------------------------------------------------------------------------------------------------------------------------------

Também são necessárias duas dependências externas:

- MGLTools-1.5.6  (gera os termos do AutoDock VIna)
- vina4dv (calcula o escore do AutoDock Vina)

É necessário criar algumas variáveis do sistema para configurar os programas externos, portanto, incluir dentro do .bashrc:


	# Set MGLTools-1.5.6
	export PATH=$PATH:/home/myuser/MGLTools-1.5.6/bin
	export MGL=/home/myuser/MGLTools-1.5.6/
	export MGLPY=$MGL/bin/python
	export MGLUTIL=$MGL/MGLToolsPckgs/AutoDockTools/Utilities24/

	# Set vina4v
	export VINADIR=/home/myuser/vina4dv/build/linux/release/

	# Set SF723_Descriptors folder
	export PATH=/home/myuser/SF723_Descriptors:$PATH



Pattern para executar SF723_Descriptors:


sf723_descriptors.py -r <path da proteína> -l <path do ligante> -f <path do arquivo FASTA> -x <lista de conjuntos de atributos> -n <nome do arquivo de saída> -o <path do arquivo de saída>


Pode gerar 8 grupos de atributos e o escore do AutoDock Vina:

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
| ID 		| NOME 							       |  Nº |
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
| amino20 	| Percentagens de aminoácidos 				       |  20 |
--------------------------------------------------------------------------------------
| dssp34 	| Descritores da estrutura secundária das proteínas 	       |  34 |
-------------------------------------------------------------------------------------
| nn350 	| Descritores do NNscore 2.0 				       | 350 |
--------------------------------------------------------------------------------------
| pd92 		| Descritores do PaDEL-Descriptors (pode ser aumentado)        |  92 |
-------------------------------------------------------------------------------------
| rdkt2d147 	| Descritores 2D gerados com RDKit 			       | 147 |
--------------------------------------------------------------------------------------
| rdkt3d11 	| Descritores 3D gerados com RDKit 			       |  11 |
--------------------------------------------------------------------------------------
| sasa10 	| Descritores da área de superfície acessível pelo solvente    |  10 |
--------------------------------------------------------------------------------------
| vina58 	| Termos do AutoDock Vina 				       |  58 |
-------------------------------------------------------------------------------------
| deltavina20 	| Escore do AutoDock Vina 				       |  20 |
--------------------------------------------------------------------------------------


São utilizados três arquivos de entrada: arquivo pdb da proteína, arquivo mol2 do ligante e arquivo FASTA da proteína com formato do The Protein Data Bank. O arquivo FASTA apenas é obligatorio para calcular a porcentagem de aminoácidos da proteína (até agora não encontrei nenhum bom conversor que calcule esse dados a partir do arquivo pdb), para os outros grupos de atributos, não é necessário.  Cada conjunto de características é criado em arquivos individuais, além de um arquivo com todos eles. Segue em anexo um vídeo de exemplos.


Alguns scripts de exemplo:

# Gerar todos os descritores

sf723_descriptors.py -r examples/1a30/1a30_protein.pdb -f examples/1a30/1A30.fasta.txt -l examples/1a30/1a30_ligand.mol2 -n 1a30 -o examples/1a30/results/ -x "amino20 dssp34 nn350 pd92 rdkt2d147 rdkt3d11 sasa10 vina58 deltavina20"

# Gerar apenas os descritores de Vina
sf723_descriptors.py -r examples/1g2k/1g2k_protein.pdb -l examples/1g2k/1g2k_ligand.mol2 -n 1g2k -o examples/1g2k/results/ -x "vina58"

# Gerar descritores de NNscore e PaDEL Descriptors
sf723_descriptors.py -r examples/1k1i/1k1i_protein.pdb -l examples/1k1i/1k1i_ligand.mol2 -n 1k1i -o examples/1k1i/results/ -x "nn350 pd92"

