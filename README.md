# Docking Molecular: da estrutura à interação

Material didático e prática em servidores web sobre docking molecular: prever onde e como uma molécula pequena se encaixa numa proteína. Dois sistemas opostos, resolvidos no mesmo servidor: o ensitrelvir na protease principal do SARS-CoV-2, em que o programa acerta a pose do cristal, e o cap do mRNA na proteína eIF4E, em que erra. Todos os arquivos de entrada estão neste repositório e baixam com um clique.

| Prática | Ferramenta | Conta | Conteúdo |
|---------|------------|-------|-------------------|
| A. Redocking do ensitrelvir | [SwissDock](https://www.swissdock.ch/) (AutoDock Vina) | não | docking, RMSD contra o cristal, interações |
| B. O cap na eIF4E | SwissDock (Vina e Attracting Cavities) | não | o mesmo protocolo falhando, e por quê |
| C. Cavidades | [PrankWeb](https://prankweb.cz/) (P2Rank) | não | achar cavidades e decidir quais são sítios |
| Interações | [PLIP](https://plip-tool.biotec.tu-dresden.de/) | não | ligações de hidrogênio, pontes salinas, empilhamento |
| Visualização | [Mol*](https://molstar.org/) no próprio PDB | não | ver a pose e o sítio em 3D |

| Autoria | |
|---|---|
| Autor | **Madson A. de Luna-Aragão** |
| Vínculos | Programa Interunidades de Pós-Graduação em Bioinformática, Universidade Federal de Minas Gerais (UFMG), Belo Horizonte, Brasil<br>Centre National de la Recherche Scientifique (CNRS), Institut de Biologie et Chimie des Protéines (IBCP), Lyon, França |
| E-mail | [madsondeluna@gmail.com](mailto:madsondeluna@gmail.com) |
| Site do curso | [madsondeluna.com/curso-docking](https://madsondeluna.com/curso-docking/) |
| Repositório | [madsondeluna/curso-docking](https://github.com/madsondeluna/curso-docking) |
| Curso de dinâmica molecular | [madsondeluna.com/curso-md](https://madsondeluna.com/curso-md/) |


## Roteiro da aula

| Bloco | Conteúdo | Material |
|-------|----------|----------|
| 1 | O que o docking responde e o que não responde | seção 1 |
| 2 | Chave e fechadura, encaixe induzido, seleção conformacional | seção 2 |
| 3 | As forças que seguram o ligante | seção 3 |
| 4 | Tipos de sítio; nem toda cavidade é um sítio | seção 4 |
| 5 | Tipos de docking: flexibilidade, busca, parceiros | seção 5 |
| 6 | Níveis de modelagem e quanto cada um custa | seção 6 |
| 7 | De onde vêm receptores e ligantes; ler uma SMILES | seção 7 |
| 8 | Métricas: energia, RMSD, interações, eficiência | seção 8 |
| 9 | Mutações, resistência e câncer | seção 9 |
| 10 | Por que o docking precisa da dinâmica molecular | seção 10 |
| 11 | Prática A, B e C nos servidores | seções de prática |


## 1. O que o docking molecular faz

O docking responde a duas perguntas sobre um par receptor e ligante: **onde** a molécula pequena se encaixa e **em que orientação e conformação**. O programa gera milhares de posições e orientações candidatas dentro de uma caixa de busca, calcula para cada uma um valor de energia com uma função de pontuação e devolve as melhores, ordenadas. Cada candidata é uma **pose**.

![Vocabulário do docking](imagens/fund_vocabulario.png)

A pergunta que o docking **não** responde bem é quanto a molécula se liga de fato. A função de pontuação é uma aproximação rápida, calibrada em complexos conhecidos; ela acerta a geometria com muito mais frequência do que a afinidade. Prever a pose e prever a afinidade são problemas diferentes, e o segundo continua difícil (Kitchen *et al.*, 2004).

**Analogia: estacionar na vaga certa.** Achar a vaga e entrar nela com o carro na posição certa é a pose. Dizer quanto tempo o carro vai ficar estacionado é a afinidade. O docking faz a primeira tarefa bem e a segunda mal.

**Exemplo saudável.** A ubiquitina, a eIF4E e milhares de outras proteínas funcionam porque reconhecem o parceiro certo entre milhares de moléculas na célula. A eIF4E segura o cap do mRNA, uma guanina metilada, e assim decide quais mensagens serão traduzidas.

**Exemplo patológico e translacional.** O ensitrelvir, antiviral oral contra a COVID-19, foi descoberto por triagem virtual seguida de desenho baseado em estrutura: um programa de docking testou uma biblioteca interna contra a protease principal do vírus, e os compostos escolhidos foram otimizados a partir das estruturas (Unoh *et al.*, 2022). É o fármaco que a Prática A usa.


## 2. Como o ligante reconhece o receptor

Três modelos explicam o encaixe, em ordem histórica, e todos continuam úteis.

![Três modelos de reconhecimento](imagens/fund_modelos_ligacao.png)

| Modelo | Ideia | Consequência para o docking |
|--------|-------|------------------------------|
| Chave e fechadura (Fischer, 1894) | o bolsão já tem a forma do ligante | basta testar posições e orientações do ligante rígido |
| Encaixe induzido (Koshland, 1958) | o bolsão muda ao receber o ligante | um receptor rígido pode não ter a forma certa |
| Seleção conformacional (Boehr, Nussinov e Wright, 2009) | o receptor oscila entre formas; o ligante escolhe uma | a forma certa pode não ser a do cristal que você baixou |

**Analogia: a luva.** Uma luva de borracha parada tem uma forma; a mão entra e ela assume outra. Se você fotografar a luva vazia e tentar encaixar a mão naquela foto rígida, pode concluir que a mão não cabe.

**Por que isso importa na prática.** Quase todo programa de docking de uso geral, inclusive o AutoDock Vina desta aula, mantém o receptor rígido. Isso funciona quando a estrutura usada já veio com um ligante parecido, e falha quando o sítio precisa se abrir. A Prática B é um caso em que isso aparece.


## 3. As forças que seguram o ligante

A energia de um complexo vem de interações fracas somadas, mais o custo de tirar a água do caminho.

![Forças que seguram o ligante](imagens/fund_forcas.png)

| Força | Alcance típico | Onde aparece nesta aula |
|-------|----------------|--------------------------|
| Ligação de hidrogênio | 2,7 a 3,2 Å entre os átomos pesados | ensitrelvir com Gly143, Ser144, Cys145 e Glu166 da Mpro |
| Interação iônica (ponte salina) | alguns ångströms, depende do meio | os três fosfatos do cap com Arg157, Lys162 e Glu103 da eIF4E |
| Efeito hidrofóbico | contato entre superfícies apolares | o ensitrelvir com Met165 |
| Empilhamento π | 3,4 a 4 Å entre anéis | o ensitrelvir com His41 |
| Cátion-π | anel aromático e carga positiva | a guanina metilada, com carga positiva, entre Trp56 e Trp102 (Ma e Dougherty, 1997) |
| Dessolvatação | custo de remover a água | penaliza ligantes muito carregados, como o cap |

**Analogia: velcro.** Nenhum gancho segura sozinho; o conjunto segura. Tirar um gancho enfraquece pouco, tirar metade solta tudo.

**Exemplo saudável.** O reconhecimento do cap pela eIF4E é um sanduíche cátion-π: a guanina metilada, carregada positivamente, fica entre dois triptofanos (Tomoo *et al.*, 2002). Essa é a interação que a Prática B tenta reproduzir.

**Exemplo patológico.** Caudas e cargas erradas quebram o reconhecimento. Na anemia falciforme, trocar um glutamato por uma valina põe um ponto hidrofóbico na superfície da hemoglobina, e as moléculas passam a se grudar umas nas outras em vez de ficarem dispersas.


## 4. Onde o ligante se encaixa: tipos de sítio

Uma cavidade é um buraco na superfície da proteína. Um sítio de ligação é uma cavidade onde algo se liga de fato. **Toda proteína tem muito mais cavidades do que sítios.**

![Tipos de sítio](imagens/fund_tipos_sitio.png)

| Tipo de sítio | O que é | Exemplo verificado nesta aula |
|---------------|---------|-------------------------------|
| Ortostérico | onde o substrato ou ligante natural se liga | sítio do ATP da quinase ABL, com nilotinibe (PDB 5MO4) |
| Alostérico | outro lugar da mesma proteína; ligar ali muda a atividade do sítio ortostérico | bolso do miristoil da ABL, com asciminibe (mesma estrutura) |
| Críptico | só se abre quando o ligante chega; não aparece, ou aparece raso, na estrutura livre | bolso switch II do KRAS G12C, com sotorasibe (PDB 6OIM) |
| Interface proteína-proteína | superfície larga e rasa que recebe outra proteína | sulco do MDM2 que recebe a hélice da p53 (PDB 1YCR) |

**Nem toda cavidade é um sítio.** Rodei o P2Rank, pelo PrankWeb, nas quatro estruturas. Ele devolve as cavidades ordenadas por um escore. Comparando o centro de cada cavidade com a posição real dos ligantes do cristal:

![Cavidades previstas e sítios reais](imagens/fund_cavidades_ranking.png)

Na Mpro são oito cavidades e só a primeira é o sítio ativo. As outras sete ficam de 8 a 26 Å do ensitrelvir: são reentrâncias da superfície sem ligante conhecido.

![Cavidades da Mpro](imagens/fund_cavidades.png)

**Sítio alostérico na quinase ABL.** Na estrutura 5MO4 a ABL tem dois fármacos ligados ao mesmo tempo. O P2Rank acerta os dois: a cavidade 1 (escore 27,9) é o sítio do ATP, a 2,1 Å do nilotinibe, e a cavidade 3 (escore 5,0) é o bolso do miristoil, a 0,9 Å do asciminibe. Um sítio alostérico é uma cavidade de verdade, só costuma ser menor e pior pontuada que a ortostérica.

**Consequência clínica de um sítio alostérico.** O asciminibe se liga ao bolso onde normalmente entra a cauda miristoilada da própria proteína e trava a quinase numa forma inativa (Wylie *et al.*, 2017). Como não disputa o sítio do ATP, mantém atividade em leucemias com a mutação T315I, que confere resistência ao imatinibe (Gorre *et al.*, 2001). Mudar de cavidade foi a saída para a resistência.

**Sítio críptico no KRAS.** O KRAS foi considerado não medicável por décadas: a superfície é lisa, sem bolsão evidente ao lado do nucleotídeo. Ostrem *et al.* (2013) mostraram que inibidores covalentes da cisteína 12 abrem um bolso novo, o switch II, e disso veio o sotorasibe (Canon *et al.*, 2019). Na estrutura com o fármaco (6OIM), o P2Rank põe esse bolso em segundo lugar. Na estrutura sem fármaco (4OBE), ele ainda aparece, mas em terceiro lugar e com escore 9,5 contra 21,5 do sítio do GDP: existe, porém raso. Uma cavidade pode ser pequena demais hoje e virar sítio amanhã.

**Cuidados ao olhar uma cavidade.**

- Buracos onde ficam glicerol, PEG, sulfato ou tampão no cristal costumam ser artefatos da cristalização, não sítios.
- Cavidades em interfaces entre cópias da proteína no cristal podem não existir na célula.
- Uma interface proteína-proteína raramente é uma cavidade funda, e por isso escapa de quem só procura bolsões.
- Um escore alto diz que a cavidade **parece** um sítio, não que algum ligante se ligue ali.


## 5. Tipos de docking

**Pela flexibilidade.**

![Tipos de docking pela flexibilidade](imagens/fund_tipos_flexibilidade.png)

| Tipo | O que se move | Custo | Programas |
|------|---------------|-------|-----------|
| Rígido-rígido | nada | muito baixo | buscas por forma, triagem inicial |
| Ligante flexível, receptor rígido | ligações rotáveis do ligante | baixo | AutoDock Vina, SwissDock, DockThor |
| Flexível-flexível | ligante e cadeias laterais, ou mais | alto | Rosetta, HADDOCK |

**Pelo que se sabe do sítio.**

![Tipos de docking pela busca](imagens/fund_tipos_busca.png)

| Tipo | Quando usar | Risco |
|------|-------------|-------|
| Baseado em molde | existe complexo parecido resolvido | depende de o molde ser mesmo parecido |
| Local | o sítio é conhecido (esta prática) | se a caixa estiver errada, a resposta está errada |
| Cego (global) | o sítio é desconhecido | muitas cavidades candidatas, mais falsos positivos |

**Pelos parceiros.** Proteína e molécula pequena é o caso desta prática. Proteína e peptídeo e proteína e proteína pedem outros programas, porque a superfície de contato é larga e rasa.

![Docking proteína-proteína](imagens/fund_barnase_barstar.png)

A barnase e a barstar formam um dos complexos mais fortes já medidos entre duas proteínas, com uma interface de dezenas de resíduos (Buckle, Schreiber e Fersht, 1994). Para esse tipo de problema existem servidores próprios, como o [ClusPro](https://cluspro.bu.edu/), que também tem modo para peptídeos (Kozakov *et al.*, 2017; Porter *et al.*, 2017).


## 6. Níveis de modelagem e quanto cada um custa

![Níveis de modelagem](imagens/fund_niveis_custo.png)

| Nível | O que modela | O que você ganha | Ordem de grandeza por ligante |
|-------|--------------|------------------|-------------------------------|
| Modelo de aprendizado de máquina | aprende poses a partir do PDB | velocidade | segundos em GPU |
| Docking com receptor rígido | forma e interações, ligante flexível | pose, ranking grosseiro | 2 min medidos nesta aula |
| Docking com cadeias laterais flexíveis | parte do encaixe induzido | sítios que precisam abrir | dezenas de minutos |
| Ensemble docking | várias conformações do receptor, vindas da dinâmica | receptor que se move | horas (Amaro *et al.*, 2018) |
| Dinâmica molecular com MM/GBSA | o complexo em água, com movimento | estabilidade da pose e energia melhor | horas a dias (Genheden e Ryde, 2015) |
| Perturbação de energia livre (FEP) | diferença de energia livre entre dois ligantes | a melhor previsão de afinidade disponível | dias em GPU (Wang *et al.*, 2015) |

**Como escolher o nível.** Use o mais barato que responda à sua pergunta. Para escolher 50 compostos entre um milhão, docking rígido basta. Para decidir qual de dois análogos sintetizar, docking não basta.

**Triagem em escala de milhões.** Lyu *et al.* (2019) testaram por docking mais de 170 milhões de moléculas que nem existiam em estoque, e compraram e testaram os melhores em bancada. Nessa escala, mesmo alguns segundos por molécula viram semanas de computador.

**Atenção com o aprendizado de máquina.** Modelos como o DiffDock produzem poses em segundos, mas uma avaliação independente mostrou que boa parte das poses geradas por métodos de IA viola geometria básica, como comprimentos de ligação e choques com a proteína (Buttenschoen, Morris e Deane, 2024). Velocidade não substitui conferir a pose.


## 7. De onde vêm receptores e ligantes

**Receptores.** O Protein Data Bank guarda estruturas resolvidas por cristalografia de raios X, ressonância magnética nuclear e crio-microscopia eletrônica. Toda estrutura tem limitações que afetam o docking:

| O que olhar na estrutura | Por que importa |
|--------------------------|-----------------|
| Resolução | abaixo de 2,5 Å as cadeias laterais do sítio são confiáveis |
| Ligante presente | uma estrutura com ligante parecido já tem o sítio aberto |
| Resíduos faltando | alças não resolvidas podem fazer parte do sítio |
| Hidrogênios | quase nunca estão no cristal; o programa os adiciona |
| Moléculas do tampão | glicerol, PEG e sulfato ocupam cavidades sem serem ligantes |

Modelos previstos, como os do AlphaFold (Jumper *et al.*, 2021), ampliaram muito a oferta de receptores, mas trazem o problema da seção 2: a conformação prevista pode não ser a que recebe o ligante.

**Ligantes.** PDB, PubChem (Kim *et al.*, 2022), ChEMBL e DrugBank fornecem estruturas e identificadores. A forma mais portátil de escrever uma molécula pequena é a SMILES, uma linha de texto (Weininger, 1988).

![A mesma molécula em desenho e em SMILES](imagens/fund_smiles.png)

![Como ler uma SMILES](imagens/fund_smiles_leitura.png)

**Como obter a SMILES de um ligante, na prática.**

1. No RCSB PDB: abra a página do ligante, por exemplo [o ligante 7YY](https://www.rcsb.org/ligand/7YY), e copie a linha SMILES isomérica.
2. No PubChem: busque pelo nome do composto e copie a SMILES canônica ou isomérica.
3. No próprio SwissDock: desenhe a molécula no editor e deixe o servidor gerar a SMILES.
4. A partir de um arquivo que você já tem: o Open Babel converte SDF, MOL2 ou PDB em SMILES.

**Protonação importa.** A mesma molécula tem cargas diferentes em pH diferente. O cap m7GTP, em pH 7,4, fica com carga total −3, por causa dos três fosfatos. Ao colar uma SMILES no SwissDock, o servidor lista os estados possíveis e pede que você escolha; o material desta aula usa sempre o estado de pH 7,4, já gravado em `dados/`.


## 8. Métricas: como ler um resultado de docking

Quatro números, e o que conta como bom e como ruim em cada um.

![Métricas do docking](imagens/fund_metricas.png)

| Métrica | O que mede | Bom | Ruim | Armadilha |
|---------|-----------|-----|------|-----------|
| Energia (kcal/mol) | pontuação da função do programa | abaixo de −9 | acima de −6 | não é energia livre medida; não compara entre programas |
| RMSD (Å) | distância média entre a pose e uma referência | abaixo de 2 | acima de 3 | só existe quando há pose de referência |
| Interações recuperadas | contatos do cristal que a pose repete | acima de 70% | abaixo de 30% | contar contato sem olhar quais contatos engana |
| Eficiência do ligante | energia dividida pelo número de átomos pesados | acima de 0,3 | abaixo de 0,2 | molécula grande ganha energia por tamanho, não por encaixe |

**Como o RMSD é calculado.** Para cada átomo pesado, mede-se a distância entre a posição na pose e a posição no cristal; o RMSD é a raiz da média dos quadrados dessas distâncias. O limiar de 2 Å é convenção da área: abaixo dele, a pose reproduz os mesmos contatos.

![Como o RMSD é calculado](imagens/fund_rmsd.png)

**Energia e constante de dissociação.** Dá para converter a pontuação numa estimativa grosseira de Kd por Kd ≈ exp(ΔG/RT), com RT ≈ 0,593 kcal/mol a 25 °C. A pose 1 da Prática A, com −9,1 kcal/mol, daria cerca de 200 nM; a da Prática B, com −6,5, daria cerca de 17 µM. Trate esses valores como ordem de grandeza, nunca como medida.

**Eficiência do ligante.** Dividir a energia pelo número de átomos pesados evita premiar moléculas grandes só por serem grandes (Hopkins, Groom e Alex, 2004). O ensitrelvir, com 37 átomos pesados, fica em 0,25 kcal/mol por átomo; o cap, com 33, fica em 0,20.

**Interações.** O PLIP lista as interações de um complexo por resíduo e por tipo (Adasme *et al.*, 2021). Comparando a pose com o cristal, dá para ver não só quantas interações se repetem, mas quais.

![Impressão digital das interações](imagens/pratica_interacoes.png)

**O gráfico que resume tudo.** Pondo energia contra RMSD para todas as poses das duas práticas, a lição aparece de uma vez: na Mpro a melhor energia é também a pose certa; na eIF4E, nenhuma pose chega perto do cristal, e a menos errada é uma das piores em energia.

![Energia contra RMSD](imagens/pratica_energia_rmsd.png)


## 9. Mutações, resistência e câncer

Uma troca de um único resíduo no sítio muda o encaixe e pode anular um fármaco.

![Imatinibe e a Thr315 da ABL](imagens/fund_resistencia_abl.png)

**Leucemia mieloide crônica.** O imatinibe se liga ao sítio do ATP da quinase BCR-ABL. A troca da treonina 315 por isoleucina, T315I, remove uma hidroxila que fazia ligação de hidrogênio com o fármaco e põe uma cadeia lateral maior no caminho. Pacientes com essa mutação deixam de responder (Gorre *et al.*, 2001; Shah *et al.*, 2002). A saída veio em duas frentes: o ponatinibe, desenhado para contornar o volume extra (O'Hare *et al.*, 2009), e o asciminibe, que troca de cavidade e vai para o sítio alostérico da seção 4.

**Câncer de pulmão.** No EGFR, a mutação T790M faz o mesmo contra o gefitinibe e o erlotinibe (Kobayashi *et al.*, 2005). É a mesma lógica: um resíduo maior no fundo do bolsão.

**COVID-19.** Sob pressão do nirmatrelvir, a Mpro do SARS-CoV-2 acumula mutações que reduzem a ligação do fármaco sem matar a enzima (Iketani *et al.*, 2022; Duan *et al.*, 2023). A estrutura 8HBK desta aula vem justamente de um estudo sobre resistência.

**Proteína-proteína e fármaco no mesmo sítio.** O MDM2 segura a p53 por um sulco que recebe três resíduos da hélice: Phe19, Trp23 e Leu26 (Kussie *et al.*, 1996). As nutlinas ocupam esse mesmo sulco e soltam a p53, reativando a resposta ao dano no DNA (Vassilev *et al.*, 2004). Um sítio que existe para receber uma proteína pode receber um fármaco.

![MDM2 com a p53 e com a nutlina](imagens/fund_p53_mdm2.png)

**O que o docking consegue dizer sobre uma mutação.** Com cuidado: comparar a pose e a energia no selvagem e no mutante dá uma indicação, mas o receptor rígido não reproduz o rearranjo que costuma acompanhar a mutação. Para isso, a seção 10.


## 10. Por que o docking precisa da dinâmica molecular

![Docking e dinâmica molecular](imagens/fund_md_docking.png)

| Momento | O que a dinâmica faz | Por quê |
|---------|---------------------|---------|
| Antes | gera várias conformações do receptor para rodar o docking em cada uma | o receptor rígido é a maior limitação do método (Lin *et al.*, 2002; Amaro *et al.*, 2018) |
| Antes | revela bolsos crípticos que não aparecem no cristal | foi assim que se encontraram mais de 50 bolsos nas proteínas do SARS-CoV-2 (Zimmerman *et al.*, 2021) |
| Depois | testa se a pose se mantém em água, a 310 K | uma pose boa no escore pode se desfazer em poucos nanossegundos |
| Depois | estima a energia com MM/GBSA ou FEP | funções de pontuação de docking são rápidas, não exatas |

**Como testar a estabilidade de uma pose.** Simule o complexo e acompanhe o RMSD do ligante ao longo do tempo. Se ele oscila numa faixa estreita, a pose é estável. Se cresce e o ligante sai do sítio, o docking entregou uma pose plausível na geometria e errada na física.

O curso de dinâmica molecular deste mesmo conjunto cobre essa parte: [madsondeluna.com/curso-md](https://madsondeluna.com/curso-md/).


## Prática A: o ensitrelvir na Mpro do SARS-CoV-2

Passo a passo com todos os cliques e valores: [guias/mpro.md](guias/mpro.md).

**Baixe as entradas:** [dados/mpro.zip](dados/mpro.zip) (receptor, ligante do cristal, SMILES em pH 7,4 e a caixa de busca).

![Os oito passos da prática](imagens/pratica_fluxo.png)

**A proteína e o fármaco.** A protease principal (Mpro, ou 3CL) corta a poliproteína do vírus em pontos específicos; sem ela, o vírus não se replica. O ensitrelvir é um inibidor não covalente, de uso oral, aprovado no Japão. A estrutura 8HBK traz os dois juntos, a 1,8 Å.

**O procedimento, em cinco passos.** 1. Conferir a estrutura no PDB. 2. Remover o ligante e guardar a posição original como gabarito. 3. Montar a caixa de busca sobre o sítio. 4. Rodar o docking só com a proteína. 5. Comparar a pose obtida com a guardada.

Esse procedimento se chama **redocking** e mede se o protocolo recupera uma resposta conhecida. É o requisito mínimo antes de aplicar o mesmo protocolo a moléculas cuja pose ninguém conhece.

| Parâmetro | Valor |
|-----------|-------|
| Receptor | 8HBK, cadeia A, 2 316 átomos pesados |
| Ligante | ensitrelvir (código 7YY), 37 átomos pesados, 5 ligações rotáveis, carga 0 |
| Centro da caixa | −21,75  −11,26  2,66 |
| Tamanho da caixa | 28 × 27 × 26 Å |
| Exaustividade | 8 |

O trabalho levou 2 minutos no SwissDock e devolveu 20 poses.

![Poses da Prática A](imagens/pratica_mpro_poses.png)

| Pose | Energia (kcal/mol) | RMSD (Å) | Interações do cristal recuperadas |
|------|--------------------|----------|-----------------------------------|
| 1 | −9,13 | 1,51 | 75% |
| 2 | −8,71 | 8,32 | 25% |
| 3 | −8,64 | 6,73 | 50% |

A pose 1 fica abaixo de 2 Å: acertou. E repete as ligações de hidrogênio com Gly143, Ser144, Cys145 e Glu166, o contato hidrofóbico com Met165 e o empilhamento com His41. A diferença de energia para a pose 2 é de apenas 0,4 kcal/mol, e a pose 2 está 8 Å fora: **a ordem do ranking é frágil, e olhar só o primeiro número é arriscado.**

**Resultados de referência.** Se o servidor estiver fora do ar no dia da aula, os arquivos que eu obtive estão em [resultados/mpro](resultados/mpro): poses, energias, RMSD e interações.


## Prática B: o cap do mRNA na eIF4E, quando o docking erra

Passo a passo: [guias/eif4e.md](guias/eif4e.md).

**Baixe as entradas:** [dados/eif4e.zip](dados/eif4e.zip).

**O sistema.** Toda molécula de mRNA eucariótico começa com um cap: uma guanina metilada ligada ao resto da fita por três fosfatos. A eIF4E reconhece esse cap e recruta a maquinaria de tradução. É o passo limitante da síntese de proteínas, e por isso a eIF4E é alvo em câncer: células tumorais com eIF4E ativa demais traduzem preferencialmente mensagens ligadas a proliferação.

**O que o cristal mostra.** Na estrutura 1IPC, a guanina metilada fica ensanduichada entre Trp56 e Trp102, numa interação cátion-π, e os fosfatos fazem pontes salinas com Arg157, Lys162 e Glu103.

Com o Vina, no mesmo protocolo da Prática A:

![Poses da Prática B](imagens/pratica_eif4e_poses.png)

| Pose | Energia (kcal/mol) | RMSD (Å) |
|------|--------------------|----------|
| 1 (melhor energia) | −6,53 | 6,51 |
| 16 (menor RMSD das 20) | −4,98 | 4,35 |

Nenhuma das 20 poses chega perto do cristal. A pose 1 mantém o empilhamento com Trp102, mas **perde o Trp56**: o sanduíche não se forma.

**Causas possíveis, e o que cada uma prevê.**

| Causa | Teste | O que aconteceu |
|-------|-------|-----------------|
| Ligante flexível demais (10 ligações rotáveis) | repetir com o m7GMP, que tem um fosfato só | continuou errando: melhor pose a 4,95 Å |
| Função de pontuação sem termo eletrostático explícito | trocar o método | melhorou, ver abaixo |
| Carga alta (−3) e dessolvatação mal tratada | o cap perde muita água ao entrar | consistente com o erro |
| Receptor rígido | o sítio da eIF4E fecha sobre o cap | não testado aqui |

**Trocando o método no mesmo servidor.** O Attracting Cavities, também do SwissDock, levou 15 minutos e devolveu 106 poses em 50 agrupamentos. A melhor pose do ranking fica em 4,60 Å, e várias poses ficam entre 2,4 e 3,0 Å. Nessas, o anel da guanina metilada cai a 0,85 Å da posição do cristal e são os fosfatos que ficam fora. O método mais caro gera a pose correta, mas a classifica em posição intermediária.

**O que a prática B estabelece.** Um protocolo que acerta num sistema pode errar em outro, e a validação por redocking vale para o alvo em que foi feita. Sem referência experimental não existe critério interno para separar uma pose correta de uma errada.


## Prática C: cavidades e sítios no PrankWeb

Passo a passo com as quatro estruturas e os critérios de julgamento: [guias/cavidades.md](guias/cavidades.md).

**O procedimento, em quatro passos.** 1. Submeter um código PDB ao [PrankWeb](https://prankweb.cz/). 2. Ler a lista de cavidades ordenadas por escore. 3. Localizar cada cavidade no visualizador e ver de qual ligante ela está perto. 4. Decidir quais cavidades são sítios e quais são reentrâncias sem função conhecida.

| Estrutura | Pergunta | Resposta medida |
|-----------|----------|------------------|
| 8HBK (Mpro) | quantas cavidades, e quantas são sítios? | 8 cavidades, 1 sítio |
| 5MO4 (ABL) | o sítio alostérico aparece? | sim, em terceiro lugar |
| 6OIM (KRAS com sotorasibe) | o bolso críptico aparece? | sim, em segundo lugar |
| 4OBE (KRAS sem fármaco) | ele aparece sem o fármaco? | aparece raso, em terceiro lugar |

Os arquivos com todas as cavidades e as distâncias até os ligantes estão em [resultados/cavidades](resultados/cavidades).

**As três partes juntas.** A Prática C trata de onde procurar, a A do método funcionando e a B do método falhando. Um resultado de docking exige as três perguntas respondidas: a cavidade é mesmo um sítio, o protocolo acerta num caso de referência, e a pose sobrevive a um teste independente.


## Referência técnica

### Estrutura do repositório

| Caminho | Conteúdo |
|---------|----------|
| `dados/<sistema>/` | receptor limpo, ligante do cristal, ligante 3D, SMILES em pH 7,4 e caixa |
| `dados/<sistema>.zip` | os mesmos arquivos num pacote, para baixar com um clique |
| `resultados/<sistema>/` | poses, energias, RMSD e interações que eu obtive nos servidores |
| `resultados/cavidades/` | cavidades previstas pelo P2Rank e distância até cada ligante |
| `guias/` | passo a passo de cada prática |
| `tools/` | scripts que preparam as entradas, rodam os servidores e geram as figuras |

### Scripts

| Script | O que faz |
|--------|-----------|
| `tools/prepare_inputs.py` | monta `dados/` a partir dos códigos PDB |
| `tools/swissdock_cli.py` | roda o docking no SwissDock pela API |
| `tools/analyze_poses.py` | energia e RMSD de cada pose contra o cristal |
| `tools/interactions.py` | interações do complexo com o PLIP |
| `tools/pockets.py` | cavidades no PrankWeb e comparação com os sítios reais |
| `tools/make_figures.py` | gera todas as figuras |
| `tools/export_results.py` | publica os resultados de referência e empacota as entradas |

### Erros comuns e o que fazer

| Sintoma | Causa provável | Solução |
|---------|----------------|---------|
| O servidor recusa o trabalho | fila cheia ou caixa grande demais | reduzir a caixa ou a exaustividade e repetir |
| Todas as poses ficam fora do sítio | caixa no lugar errado | conferir o centro contra o ligante do cristal |
| Energias muito melhores que o esperado | ligante grande demais ou com clash | olhar a pose, conferir a eficiência do ligante |
| O ligante sai distorcido | SMILES com estado de protonação errado | usar o estado de pH 7,4 |
| RMSD impossível de calcular | falta a pose de referência | só dá para calcular em redocking |

### Ferramentas usadas

| Ferramenta | Papel | Acesso |
|------------|-------|--------|
| SwissDock (Vina e Attracting Cavities) | docking | livre, sem conta |
| PrankWeb (P2Rank) | cavidades | livre, sem conta |
| PLIP | interações | livre, sem conta |
| Mol* | visualização | livre, sem conta |
| ClusPro | proteína-proteína e peptídeo | livre, uso não comercial, conta opcional |
| DiffDock (NVIDIA) | docking por aprendizado de máquina | exige conta e chave de API; não testado neste material |


## Referências

- Trott, O.; Olson, A. J. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *J. Comput. Chem.* **2010**, *31*, 455-461. https://doi.org/10.1002/jcc.21334
- Bugnon, M. *et al.* SwissDock 2024: major enhancements for small-molecule docking with Attracting Cavities and AutoDock Vina. *Nucleic Acids Res.* **2024**, *52*, W324-W332. https://doi.org/10.1093/nar/gkae300
- Jakubec, D. *et al.* PrankWeb 3: accelerated ligand-binding site predictions for experimental and modelled protein structures. *Nucleic Acids Res.* **2022**, *50*, W593-W597. https://doi.org/10.1093/nar/gkac389
- Jendele, L. *et al.* PrankWeb: a web server for ligand binding site prediction and visualization. *Nucleic Acids Res.* **2019**, *47*, W345-W349. https://doi.org/10.1093/nar/gkz424
- Le Guilloux, V.; Schmidtke, P.; Tuffery, P. Fpocket: an open source platform for ligand pocket detection. *BMC Bioinformatics* **2009**, *10*, 168. https://doi.org/10.1186/1471-2105-10-168
- Adasme, M. F. *et al.* PLIP 2021: expanding the scope of the protein-ligand interaction profiler to DNA and RNA. *Nucleic Acids Res.* **2021**, *49*, W530-W534. https://doi.org/10.1093/nar/gkab294
- Sehnal, D. *et al.* Mol* Viewer: modern web app for 3D visualization and analysis of large biomolecular structures. *Nucleic Acids Res.* **2021**, *49*, W431-W437. https://doi.org/10.1093/nar/gkab314
- Kozakov, D. *et al.* The ClusPro web server for protein-protein docking. *Nat. Protoc.* **2017**, *12*, 255-278. https://doi.org/10.1038/nprot.2016.169
- Porter, K. A. *et al.* ClusPro PeptiDock: efficient global docking of peptide recognition motifs using FFT. *Bioinformatics* **2017**, *33*, 3299-3301. https://doi.org/10.1093/bioinformatics/btx216
- Dominguez, C.; Boelens, R.; Bonvin, A. M. J. J. HADDOCK: a protein-protein docking approach based on biochemical or biophysical information. *J. Am. Chem. Soc.* **2003**, *125*, 1731-1737. https://doi.org/10.1021/ja026939x
- Kitchen, D. B.; Decornez, H.; Furr, J. R.; Bajorath, J. Docking and scoring in virtual screening for drug discovery: methods and applications. *Nat. Rev. Drug Discov.* **2004**, *3*, 935-949. https://doi.org/10.1038/nrd1549
- Buttenschoen, M.; Morris, G. M.; Deane, C. M. PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences. *Chem. Sci.* **2024**, *15*, 3130-3139. https://doi.org/10.1039/D3SC04185A
- Lyu, J. *et al.* Ultra-large library docking for discovering new chemotypes. *Nature* **2019**, *566*, 224-229. https://doi.org/10.1038/s41586-019-0917-9
- Fischer, E. Einfluss der Configuration auf die Wirkung der Enzyme. *Ber. Dtsch. Chem. Ges.* **1894**, *27*, 2985-2993. https://doi.org/10.1002/cber.18940270364
- Koshland, D. E. Application of a theory of enzyme specificity to protein synthesis. *Proc. Natl. Acad. Sci. USA* **1958**, *44*, 98-104. https://doi.org/10.1073/pnas.44.2.98
- Boehr, D. D.; Nussinov, R.; Wright, P. E. The role of dynamic conformational ensembles in biomolecular recognition. *Nat. Chem. Biol.* **2009**, *5*, 789-796. https://doi.org/10.1038/nchembio.232
- Ma, J. C.; Dougherty, D. A. The cation-π interaction. *Chem. Rev.* **1997**, *97*, 1303-1324. https://doi.org/10.1021/cr9603744
- Tomoo, K. *et al.* Crystal structures of 7-methylguanosine 5'-triphosphate (m7GTP)- and P1-7-methylguanosine-P3-adenosine-5',5'-triphosphate (m7GpppA)-bound human full-length eukaryotic initiation factor 4E. *Biochem. J.* **2002**, *362*, 539-544. https://doi.org/10.1042/bj3620539
- Unoh, Y. *et al.* Discovery of S-217622, a noncovalent oral SARS-CoV-2 3CL protease inhibitor clinical candidate for treating COVID-19. *J. Med. Chem.* **2022**, *65*, 6499-6512. https://doi.org/10.1021/acs.jmedchem.2c00117
- Iketani, S. *et al.* Multiple pathways for SARS-CoV-2 resistance to nirmatrelvir. *Nature* **2022**, *613*, 558-564. https://doi.org/10.1038/s41586-022-05514-2
- Wylie, A. A. *et al.* The allosteric inhibitor ABL001 enables dual targeting of BCR-ABL1. *Nature* **2017**, *543*, 733-737. https://doi.org/10.1038/nature21702
- Gorre, M. E. *et al.* Clinical resistance to STI-571 cancer therapy caused by BCR-ABL gene mutation or amplification. *Science* **2001**, *293*, 876-880. https://doi.org/10.1126/science.1062538
- Shah, N. P. *et al.* Multiple BCR-ABL kinase domain mutations confer polyclonal resistance to the tyrosine kinase inhibitor imatinib (STI571) in chronic phase and blast crisis chronic myeloid leukemia. *Cancer Cell* **2002**, *2*, 117-125. https://doi.org/10.1016/S1535-6108(02)00096-X
- O'Hare, T. *et al.* AP24534, a pan-BCR-ABL inhibitor for chronic myeloid leukemia, potently inhibits the T315I mutant and overcomes mutation-based resistance. *Cancer Cell* **2009**, *16*, 401-412. https://doi.org/10.1016/j.ccr.2009.09.028
- Kobayashi, S. *et al.* EGFR mutation and resistance of non-small-cell lung cancer to gefitinib. *N. Engl. J. Med.* **2005**, *352*, 786-792. https://doi.org/10.1056/NEJMoa044238
- Ostrem, J. M. *et al.* K-Ras(G12C) inhibitors allosterically control GTP affinity and effector interactions. *Nature* **2013**, *503*, 548-551. https://doi.org/10.1038/nature12796
- Canon, J. *et al.* The clinical KRAS(G12C) inhibitor AMG 510 drives anti-tumour immunity. *Nature* **2019**, *575*, 217-223. https://doi.org/10.1038/s41586-019-1694-1
- Kussie, P. H. *et al.* Structure of the MDM2 oncoprotein bound to the p53 tumor suppressor transactivation domain. *Science* **1996**, *274*, 948-953. https://doi.org/10.1126/science.274.5289.948
- Vassilev, L. T. *et al.* In vivo activation of the p53 pathway by small-molecule antagonists of MDM2. *Science* **2004**, *303*, 844-848. https://doi.org/10.1126/science.1092472
- Buckle, A. M.; Schreiber, G.; Fersht, A. R. Protein-protein recognition: crystal structural analysis of a barnase-barstar complex at 2.0-Å resolution. *Biochemistry* **1994**, *33*, 8878-8889. https://doi.org/10.1021/bi00196a004
- Lin, J.-H. *et al.* Computational drug design accommodating receptor flexibility: the relaxed complex scheme. *J. Am. Chem. Soc.* **2002**, *124*, 5632-5633. https://doi.org/10.1021/ja0260162
- Amaro, R. E. *et al.* Ensemble docking in drug discovery. *Biophys. J.* **2018**, *114*, 2271-2278. https://doi.org/10.1016/j.bpj.2018.02.038
- Genheden, S.; Ryde, U. The MM/PBSA and MM/GBSA methods to estimate ligand-binding affinities. *Expert Opin. Drug Discov.* **2015**, *10*, 449-461. https://doi.org/10.1517/17460441.2015.1032936
- Wang, L. *et al.* Accurate and reliable prediction of relative ligand binding potency in prospective drug discovery by way of a modern free-energy calculation protocol and force field. *J. Am. Chem. Soc.* **2015**, *137*, 2695-2703. https://doi.org/10.1021/ja512751q
- Zimmerman, M. I. *et al.* SARS-CoV-2 simulations go exascale to predict dramatic spike opening and cryptic pockets across the proteome. *Nat. Chem.* **2021**, *13*, 651-659. https://doi.org/10.1038/s41557-021-00707-0
- Jumper, J. *et al.* Highly accurate protein structure prediction with AlphaFold. *Nature* **2021**, *596*, 583-589. https://doi.org/10.1038/s41586-021-03819-2
- Kim, S. *et al.* PubChem 2023 update. *Nucleic Acids Res.* **2022**, *51*, D1373-D1380. https://doi.org/10.1093/nar/gkac956
- Weininger, D. SMILES, a chemical language and information system. 1. Introduction to methodology and encoding rules. *J. Chem. Inf. Comput. Sci.* **1988**, *28*, 31-36. https://doi.org/10.1021/ci00057a005
- Hopkins, A. L.; Groom, C. R.; Alex, A. Ligand efficiency: a useful metric for lead selection. *Drug Discov. Today* **2004**, *9*, 430-431. https://doi.org/10.1016/S1359-6446(04)03069-7


## Licença

Este material está sob a licença MIT. Veja [LICENSE](LICENSE).
