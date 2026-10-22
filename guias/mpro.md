# Passo a passo: o ensitrelvir na Mpro do SARS-CoV-2

Prática A, do começo ao fim. Cada passo diz o que abrir, o que clicar, o que digitar, por que aquela escolha e o que deve aparecer na tela. Os valores de referência vêm da execução guardada em `resultados/mpro`.

[Voltar à aula](../README.md) · [Prática B: eIF4E](eif4e.md) · [Prática C: cavidades](cavidades.md)

![Os oito passos da prática](../imagens/pratica_fluxo.png)

## O que será feito, em uma frase

A estrutura 8HBK traz a protease do SARS-CoV-2 com o fármaco ensitrelvir encaixado. Você vai **apagar o fármaco**, entregar ao programa só a proteína e pedir que ele descubra sozinho onde a molécula se encaixa. No fim, compara a resposta do programa com a posição verdadeira, que ficou guardada.

Esse teste tem nome: **redocking**. Serve para medir se o protocolo acerta num caso em que a resposta é conhecida. Só depois de passar nele é que faz sentido usar o mesmo protocolo em moléculas novas, cuja resposta ninguém sabe.

## Antes de começar

Precisa só de um navegador. Nenhum dos serviços pede cadastro.

**Baixe o pacote de entradas:** [dados/mpro.zip](../dados/mpro.zip) e descompacte.

| Arquivo | O que é | Usado no passo |
|---------|---------|----------------|
| `receptor.pdb` | cadeia A da 8HBK, só proteína | 3 |
| `ligante.smi` | SMILES do ensitrelvir em pH 7,4 | 2 |
| `caixa.json` | centro e tamanho da caixa | 4 |
| `ligante_cristal.sdf` | ensitrelvir na posição verdadeira | 7 |
| `ligante_3d.sdf` e `.mol2` | ligante em 3D com hidrogênios | alternativa ao passo 2 |

## Passo 1: inspecionar a estrutura no PDB

**1.1.** Abra [rcsb.org/structure/8HBK](https://www.rcsb.org/structure/8HBK).

**1.2.** Na caixa **Experimental Data Snapshot**, no alto à direita, confira:

| Campo | Valor esperado | Por que você confere |
|-------|----------------|----------------------|
| Method | X-RAY DIFFRACTION | as posições dos átomos são medidas, não previstas |
| Resolution | 1.60 a 1.80 Å | abaixo de 2,5 Å as cadeias laterais do sítio são confiáveis |

**1.3.** Role até **Macromolecules** e confirme que existe uma cadeia A de proteína. Em **Small Molecules**, confirme o ligante de código **7YY**: é o ensitrelvir.

**1.4.** Clique em **3D View**. Abre o visualizador Mol*. Gire com o botão esquerdo do mouse e aproxime com a roda até ver o ligante, em bastões, dentro de uma fenda da proteína.

**Por que este passo existe.** O programa de docking não sinaliza quando a estrutura é inadequada: ele devolve energias em qualquer caso. Resolução baixa, alças não resolvidas ou sítio fechado produzem resultados sem significado e com aparência normal.

**O que você deve ver.** Uma molécula alongada ocupando uma fenda entre dois domínios da proteína. Essa fenda é o sítio ativo, e é ali que a caixa do passo 4 vai ficar.

![Vocabulário do docking na própria 8HBK](../imagens/fund_vocabulario.png)

## Passo 2: informar o ligante no SwissDock

**2.1.** Abra [swissdock.ch](https://www.swissdock.ch/) e vá para a página de submissão.

**2.2.** No campo do **ligante**, entre as opções de SMILES, arquivo Mol2 e desenho, escolha **SMILES**.

**2.3.** Abra o arquivo `ligante.smi` num editor de texto e copie a linha, ou copie daqui:

```
Cn1cnc(Cn2c(=O)nc(Nc3cc4cn(C)nc4cc3Cl)n(Cc3cc(F)c(F)cc3F)c2=O)n1
```

**2.4.** Cole no campo e confirme.

**2.5.** O servidor responde com uma lista de **estados de protonação**, cada um com a faixa de pH em que predomina. **Escolha o da faixa que contém 7,4.** Para o ensitrelvir, esse estado tem carga total 0.

**Por que a protonação é uma escolha de verdade.** A mesma molécula tem cargas diferentes conforme o pH. Carga muda quais interações são possíveis: um grupo que deveria fazer ponte salina, se entrar neutro, simplesmente não faz. O pH 7,4 é o do sangue e do meio extracelular, e é o padrão quando se quer representar o ambiente fisiológico.

**Alternativa ao SMILES.** Se preferir enviar arquivo, use `ligante_3d.mol2`. O servidor exige um arquivo com **uma só molécula, com todos os hidrogênios e coordenadas 3D**, e esse arquivo já atende às três condições.

**O que esperar.** A mensagem de que o ligante foi preparado e um **número de sessão**. Anote esse número: é por ele que você acompanha e recupera o trabalho.

## Passo 3: informar o receptor

**3.1.** No campo do **alvo**, escolha enviar arquivo e selecione `receptor.pdb`.

**3.2.** Confirme o envio.

**Por que enviar o arquivo em vez de digitar o código 8HBK.** O servidor aceita o código e deixa você escolher cadeias e heteroátomos ali mesmo. Funciona, mas o arquivo do pacote já foi limpo da mesma forma para todos, o que torna os resultados comparáveis entre os alunos e evita erro de seleção.

**O que já foi retirado desse arquivo, e por quê.**

| Removido | Motivo |
|----------|--------|
| O ligante 7YY | se ficasse, não haveria espaço livre e o teste perderia o sentido |
| Águas | o programa trata o solvente de forma implícita; águas do cristal ocupariam o sítio |
| Íons, glicerol, PEG e sulfato | são aditivos da cristalização, não fazem parte da proteína |

**Se você optar pelo código PDB.** Marque só a cadeia A e **não selecione solvente químico** entre os heteroátomos. Desmarque também o próprio ligante 7YY; sem isso, não há redocking.

**O que esperar.** O visualizador do servidor passa a mostrar a proteína.

## Passo 4: definir a caixa de busca

**4.1.** Abra o arquivo `caixa.json`. Ele tem os seis números a digitar.

**4.2.** Preencha os campos:

| Campo | Valor |
|-------|-------|
| Centro X | −21,75 |
| Centro Y | −11,26 |
| Centro Z | 2,66 |
| Tamanho X | 28 |
| Tamanho Y | 27 |
| Tamanho Z | 26 |

**4.3.** Olhe o visualizador 3D ao lado. A caixa aparece desenhada sobre a proteína. **Confirme que ela envolve a fenda que você viu no passo 1.4.** Se estiver em outro lugar, confira os sinais negativos.

**De onde vieram esses números.** O centro é o centro geométrico do ligante no cristal. O tamanho é a extensão do ligante mais 8 Å de folga em cada direção, para o programa poder testar orientações sem esbarrar na parede da caixa.

**Por que a caixa é a decisão mais importante da prática.** Ela define onde o programa procura. Uma caixa fora do sítio produz poses sempre erradas, com energias de aparência normal. Uma caixa grande demais aumenta o tempo e enche o resultado de poses inúteis. Como aqui o sítio é conhecido, isto é um **docking local**.

**O que esperar.** Uma caixa que cobre o sítio ativo inteiro e sobra pouco para fora.

## Passo 5: escolher o método e enviar

**5.1.** Em amostragem, deixe a exaustividade em **8**, que é o valor da execução de referência.

**5.2.** Clique no botão **Docking with AutoDock Vina**.

**Por que o Vina e não o outro botão.** O Vina é rápido, é o programa de docking mais usado e trata o ligante como flexível e o receptor como rígido (Trott e Olson, 2010). Esse é o compromisso padrão da área: mexer no ligante é barato, mexer no receptor é caro. O outro método do servidor, Attracting Cavities, é mais caro e entra na [Prática B](eif4e.md).

**O que o servidor faz antes de aceitar.** Ele estima o tempo do trabalho e recusa se for longo demais. **Se recusar:** reduza o tamanho da caixa ou a exaustividade e envie de novo.

**O que esperar.** Uma página de estado. O trabalho de referência levou **2 minutos**.

## Passo 6: baixar e ler o resultado

**6.1.** Quando a página indicar que terminou, baixe o pacote de resultados.

**6.2.** Descompacte. Dentro dele:

| Arquivo | Conteúdo |
|---------|----------|
| `vina_dock.pdbqt` | as 20 poses com as energias |
| `vina_dock.pdb` | as mesmas poses, em formato PDB |
| `receptor.pdb` | o receptor como o servidor o usou |
| `parameters` | os valores que você preencheu, para registro |

**6.3.** Abra `vina_dock.pdbqt` num editor de texto. Cada bloco que começa em `MODEL` é uma pose, em ordem de energia. A linha logo abaixo traz a energia:

```
MODEL 1
REMARK VINA RESULT:    -9.127      0.000      0.000
```

O primeiro número é a energia em kcal/mol. Anote os valores das três primeiras poses.

**O que você deve obter.**

| Pose | Energia (kcal/mol) | RMSD contra o cristal (Å) |
|------|--------------------|----------------------------|
| 1 | −9,13 | 1,51 |
| 2 | −8,71 | 8,32 |
| 3 | −8,64 | 6,73 |
| 4 | −8,43 | 7,27 |

**Por que seus números vão diferir um pouco.** O Vina sorteia a semente do gerador aleatório, então duas execuções iguais dão resultados próximos, não idênticos. Espere a pose 1 entre −9,5 e −8,5 kcal/mol. **Se a sua pose 1 passar de 3 Å, volte ao passo 4 e confira a caixa.**

**O que fazer com a energia.** Tratar como ordem de grandeza. Pela relação Kd ≈ exp(ΔG/RT), com RT ≈ 0,593 kcal/mol, os −9,1 kcal/mol dariam cerca de 200 nM. A afinidade medida em bancada não sai daqui.

## Passo 7: comparar a pose com o cristal

**7.1.** Abra [molstar.org/viewer](https://molstar.org/viewer/).

**7.2.** No painel da esquerda, clique em **Open Files** e selecione dois arquivos: `ligante_cristal.sdf`, do pacote de entradas, e `vina_dock.pdb`, do pacote de resultados.

**7.3.** Dê cores diferentes aos dois e gire a cena. Veja se os anéis da pose caem sobre os anéis do cristal.

**O que o número significa.** O RMSD é a distância média entre cada átomo pesado da pose e o mesmo átomo no cristal. Não se realinha nada antes: as duas estruturas já estão no mesmo sistema de coordenadas, porque vieram da mesma estrutura.

![Como o RMSD é calculado](../imagens/fund_rmsd.png)

| RMSD | Leitura |
|------|---------|
| abaixo de 2 Å | pose correta, repete os mesmos contatos |
| 2 a 3 Å | orientação geral certa, detalhes errados |
| acima de 3 Å | pose errada |

**Por que o limiar é 2 Å.** É convenção da área, com razão física: acima disso os grupos funcionais deixam de apontar para os mesmos resíduos, e as interações mudam.

A pose 1 cai praticamente sobre o cristal. A pose 2, com energia apenas 0,4 kcal/mol pior, está 8 Å fora.

![Pose 1 e pose 2 contra o cristal](../imagens/pratica_mpro_poses.png)

**Por que isso só existe em redocking.** Sem uma pose experimental de referência não há RMSD para calcular. Com uma molécula nova você fica só com a energia, que é o número menos confiável.

## Passo 8: ver as interações no PLIP

**8.1. Primeiro, as interações do cristal.** Abra o [PLIP](https://plip-tool.biotec.tu-dresden.de/), digite **8HBK** no campo de código PDB e envie.

**8.2.** Na saída, procure a tabela do ligante 7YY. Você deve encontrar:

| Tipo | Resíduos |
|------|----------|
| Ligação de hidrogênio | Thr26, Gly143, Ser144, Cys145, Glu166 |
| Hidrofóbica | Thr25, Met165 |
| Empilhamento π | His41 |

**Por que esses resíduos e não outros.** Cys145 e His41 são a díade catalítica da protease: o ensitrelvir ocupa justamente o ponto onde o substrato seria cortado. Gly143, Ser144 e Cys145 formam o buraco de oxiânion, que estabiliza o intermediário da reação. O fármaco funciona porque imita as interações do substrato.

**8.3. Agora, as interações da sua pose.** O PLIP precisa de um arquivo com a proteína e o ligante juntos, e o ligante marcado como `HETATM`. Há dois caminhos.

**Caminho curto.** Baixe o arquivo pronto [resultados/mpro/complexo_pose1.pdb](../resultados/mpro/complexo_pose1.pdb) e envie ao PLIP pelo botão de upload.

**Caminho manual, com a sua própria pose.**

1. Abra `receptor.pdb` num editor de texto simples e apague a última linha, `END`.
2. Abra `vina_dock.pdb` e copie as linhas da primeira pose: começam depois de `MODEL 1` e terminam antes do primeiro `ENDMDL`. São 38 linhas que começam com `ATOM`.
3. Cole essas linhas no fim do `receptor.pdb`.
4. Nas linhas coladas, troque `ATOM  ` por `HETATM` no começo de cada uma. Sem isso o PLIP trata o ligante como parte da proteína e não analisa nada.
5. Acrescente uma linha `END` no fim e salve como `complexo.pdb`.
6. Envie esse arquivo ao PLIP.

**8.4.** Compare as duas listas. A pose 1 da referência repete 75% das interações do cristal: mantém as ligações de hidrogênio com Gly143, Ser144, Cys145 e Glu166, o contato hidrofóbico com Met165 e o empilhamento com His41. Perde só os contatos com Thr25 e Thr26, que ficam na borda do sítio.

![Impressão digital das interações](../imagens/pratica_interacoes.png)

**Por que contar interações, se já existe o RMSD.** O RMSD é um número só e trata todos os átomos por igual. Duas poses com RMSD parecido podem fazer contatos diferentes. A lista de interações diz **quais** contatos a pose acerta, e é disso que depende a atividade de um fármaco.

## O que concluir desta prática

1. A pose de melhor energia é também a pose certa: neste sistema o protocolo funciona.
2. A diferença entre a pose 1 e a pose 2 é de 0,4 kcal/mol, e a pose 2 está 8 Å fora. A ordem do ranking é frágil; olhar só o primeiro número é arriscado.
3. Um redocking bem-sucedido valida o protocolo para moléculas parecidas naquele sítio. Não valida para outro alvo, como mostra a [Prática B](eif4e.md).

## Se o servidor estiver fora do ar

Tudo o que você precisa para os passos 6, 7 e 8 está em [resultados/mpro](../resultados/mpro):

| Arquivo | Para que serve |
|---------|----------------|
| `poses.json` | energia e RMSD de cada uma das 20 poses |
| `poses_top3.sdf` | as três melhores poses em 3D, para abrir no Mol* |
| `complexo_pose1.pdb` | receptor e pose 1 juntos, pronto para o PLIP |
| `interacoes.json` | interações do cristal e das poses, já calculadas |
