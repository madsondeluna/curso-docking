# Passo a passo: o ensitrelvir na Mpro do SARS-CoV-2

Redocking do ensitrelvir na protease principal do SARS-CoV-2 (PDB 8HBK), no SwissDock, com AutoDock Vina. Cada passo traz o que fazer, o que selecionar, por que aquela escolha e o que esperar na tela. Os valores de referência vêm da execução guardada em `resultados/mpro`.

[Voltar à aula](../README.md) · [Passo a passo da eIF4E](eif4e.md) · [Passo a passo das cavidades](cavidades.md)

## Antes de começar

Só é preciso um navegador. Nenhum serviço desta prática pede cadastro.

**Baixe o pacote de entradas:** [dados/mpro.zip](../dados/mpro.zip)

| Arquivo | O que é | Quando você usa |
|---------|---------|------------------|
| `receptor.pdb` | cadeia A da 8HBK, só proteína | passo 3, campo do alvo |
| `ligante.smi` | SMILES do ensitrelvir em pH 7,4 | passo 2, campo do ligante |
| `ligante_cristal.sdf` | ensitrelvir na pose do cristal | passo 6, comparação |
| `caixa.json` | centro e tamanho da caixa | passo 4 |
| `ligante_3d.sdf` e `.mol2` | ligante em 3D com hidrogênios | alternativa ao SMILES |

**O objetivo.** Tirar o ligante que já está no cristal, rodar o docking só com a proteína e ver se o programa recoloca a molécula onde ela estava. Isso é um **redocking**, e é o teste que valida um protocolo antes de aplicá-lo a moléculas novas. Sem esse teste, não há como saber se um resultado de docking merece confiança.

## Passo 1: inspecionar a estrutura

**O que fazer.** Abra [rcsb.org/structure/8HBK](https://www.rcsb.org/structure/8HBK).

**O que olhar, nesta ordem.**

| Campo na página | Valor esperado | Por que importa |
|-----------------|----------------|------------------|
| Method | X-RAY DIFFRACTION | a posição dos átomos é experimental, não prevista |
| Resolution | 1,8 Å | abaixo de 2,5 Å as cadeias laterais do sítio são confiáveis |
| Ligands | 7YY | existe um ligante no sítio, então o bolsão está aberto |
| Organism | SARS-CoV-2 | confirma que é a proteína certa |

**O que fazer na sequência.** Clique em **3D View** para abrir o Mol* e gire a estrutura até ver o ligante dentro da fenda.

**Por que este passo existe.** O docking não avisa que a estrutura é ruim; ele devolve números de qualquer jeito. Resolução baixa, alças faltando ou sítio fechado produzem resultados sem sentido que parecem normais. Conferir a estrutura é parte do método, não formalidade.

**O que esperar.** Um ligante alongado, ocupando uma fenda entre dois domínios.

## Passo 2: informar o ligante

**O que fazer.** Abra [swissdock.ch](https://www.swissdock.ch/) e vá ao campo do ligante.

**O que selecionar.** Entre as três opções, SMILES, arquivo Mol2 ou desenho, use **SMILES** e cole a linha de `ligante.smi`:

```
Cn1cnc(Cn2c(=O)nc(Nc3cc4cn(C)nc4cc3Cl)n(Cc3cc(F)c(F)cc3F)c2=O)n1
```

**O que o servidor pergunta em seguida.** Ele lista os estados de protonação possíveis, cada um com a faixa de pH em que predomina. **Selecione o estado de pH 7,4.** No ensitrelvir, a carga total é 0.

**Por que a protonação é uma escolha, e não um detalhe.** A mesma molécula tem cargas diferentes em pH diferente, e carga muda quais interações são possíveis. Um grupo que deveria fazer ponte salina e entra neutro simplesmente não faz. Escolher o estado fisiológico, pH 7,4, é o padrão quando se quer representar o ambiente da célula.

**Alternativa.** Se preferir enviar arquivo, use `ligante_3d.mol2`. Ele já tem hidrogênios e coordenadas 3D, que é o que o servidor exige.

**O que esperar.** A mensagem de que o ligante foi preparado e um número de sessão.

## Passo 3: informar o alvo

**O que fazer.** No campo do alvo, envie o arquivo `receptor.pdb`.

**Por que não usar o código PDB direto.** O servidor aceita um código do PDB e deixa você escolher cadeias e heteroátomos. Funciona, mas o arquivo do pacote já está limpo da mesma forma para todos, o que torna os resultados comparáveis entre os alunos. Ele contém só os átomos da cadeia A, sem água, sem íons e, principalmente, **sem o ligante**.

**Por que remover a água e o ligante.**

| O que é removido | Motivo |
|------------------|--------|
| O ligante | se ficasse, não haveria espaço livre e o teste perderia o sentido |
| Águas | o programa trata o solvente de forma implícita; águas do cristal ocupariam o sítio |
| Íons e tampão | glicerol, PEG e sulfato são artefatos da cristalização, não fazem parte da proteína |
| Cadeia B, quando existe | uma cópia basta e o trabalho fica mais rápido |

**Atenção.** Se usar o código PDB em vez do arquivo, não selecione solvente químico entre os heteroátomos.

## Passo 4: definir a caixa de busca

**O que fazer.** Preencha os campos de centro e tamanho com os valores de `caixa.json`:

| Campo | Valor |
|-------|-------|
| Centro X | −21,75 |
| Centro Y | −11,26 |
| Centro Z | 2,66 |
| Tamanho X | 28 |
| Tamanho Y | 27 |
| Tamanho Z | 26 |

**Por que esses números.** O centro é o centro geométrico do ligante no cristal, e o tamanho é a extensão do ligante mais 8 Å de folga em cada direção. A folga dá espaço para o programa testar orientações e conformações sem forçar a molécula contra a parede da caixa.

**Por que a caixa é a decisão mais importante da prática.** A caixa define onde o programa procura. Uma caixa fora do sítio produz poses sempre erradas, com energias que parecem razoáveis. Uma caixa grande demais aumenta o tempo e o número de poses inúteis. Como aqui o sítio é conhecido, isto é um **docking local**: a caixa cobre só o sítio ativo.

**O que esperar.** O visualizador 3D ao lado mostra a caixa sobre a proteína. Confira que ela envolve a fenda que você viu no passo 1.

## Passo 5: escolher o método e enviar

**O que selecionar.** Clique em **Docking with AutoDock Vina**.

**Por que o Vina.** É rápido, é o programa de docking mais usado e mais citado, e trata o ligante como flexível e o receptor como rígido (Trott e Olson, 2010). Esse é o compromisso padrão da área: a flexibilidade do ligante é barata, a do receptor é cara.

**O que ajustar em amostragem.** A exaustividade controla quanto o programa procura. Mais amostragem dá resultado mais reprodutível e demora mais. A execução de referência usou exaustividade 8.

**O que o servidor faz antes de aceitar.** Ele estima o tempo e recusa trabalhos longos demais. Se recusar, reduza o tamanho da caixa ou a exaustividade e envie de novo.

**O que esperar.** Uma página de estado, e cerca de **2 minutos** até o fim. Guarde o número da sessão.

## Passo 6: ler os resultados

**O que fazer.** Baixe o pacote de resultados. O arquivo central é `vina_dock.pdbqt`.

**Como ele é organizado.** Cada bloco `MODEL` é uma pose, em ordem de energia. A linha `REMARK VINA RESULT` traz três números; o primeiro é a energia em kcal/mol.

**O que você deve obter.**

| Pose | Energia (kcal/mol) | RMSD contra o cristal (Å) |
|------|--------------------|----------------------------|
| 1 | −9,13 | 1,51 |
| 2 | −8,71 | 8,32 |
| 3 | −8,64 | 6,73 |
| 4 | −8,43 | 7,27 |

**Por que seus números vão variar um pouco.** O Vina sorteia a semente do gerador aleatório. Espere a pose 1 entre −9,5 e −8,5 kcal/mol, e abaixo de 2 Å. Se a sua pose 1 estiver acima de 3 Å, confira a caixa do passo 4.

**O que fazer com a energia.** Trate como ordem de grandeza. Pela relação Kd ≈ exp(ΔG/RT), −9,1 kcal/mol daria cerca de 200 nM. A afinidade medida em bancada não vem daqui.

## Passo 7: comparar com o cristal

**O que fazer.** Abra `ligante_cristal.sdf` e a pose 1 no mesmo visualizador, por exemplo o Mol*, e veja se os anéis caem no mesmo lugar.

**O que o número significa.** O RMSD mede a distância média entre cada átomo pesado da pose e o mesmo átomo no cristal, sem realinhar as estruturas: as duas já estão no mesmo referencial.

| RMSD | Leitura |
|------|---------|
| abaixo de 2 Å | pose correta, repete os mesmos contatos |
| 2 a 3 Å | quase, orientação geral certa e detalhes errados |
| acima de 3 Å | pose errada |

**Por que 2 Å.** É a convenção da área, e tem motivo físico: acima disso os grupos funcionais deixam de apontar para os mesmos resíduos, e as interações mudam.

**Por que isso só existe em redocking.** Sem uma pose experimental de referência, não há RMSD para calcular. Com uma molécula nova, você fica só com a energia, que é o número menos confiável.

## Passo 8: olhar as interações no PLIP

**O que fazer.** Abra o [PLIP](https://plip-tool.biotec.tu-dresden.de/) e envie o código **8HBK** para ver as interações do cristal.

**O que esperar na saída.**

| Tipo | Resíduos |
|------|----------|
| Ligação de hidrogênio | Thr26, Gly143, Ser144, Cys145, Glu166 |
| Hidrofóbica | Thr25, Met165 |
| Empilhamento π | His41 |

**Por que esses resíduos e não outros.** Cys145 e His41 formam a díade catalítica da protease: o ensitrelvir ocupa justamente o lugar onde o substrato seria cortado. Gly143, Ser144 e Cys145 formam o buraco de oxiânion, que estabiliza o intermediário da reação. O fármaco funciona porque imita as interações do substrato.

**Como comparar com a sua pose.** Monte um PDB com o receptor e a pose juntos e envie ao PLIP. A pose 1 da referência repete 75% das interações do cristal: perde apenas os contatos com Thr25 e Thr26, na borda do sítio.

**Por que contar interações, se já existe o RMSD.** O RMSD é um número só, e trata todos os átomos igualmente. Duas poses com RMSD parecido podem fazer contatos diferentes. A lista de interações diz **quais** contatos a pose acerta, e é isso que explica a atividade de um fármaco.

## O que concluir

- A pose de melhor energia é a pose certa: neste sistema o protocolo funciona.
- A diferença entre a pose 1 e a pose 2 é de 0,4 kcal/mol, e a pose 2 está 8 Å fora. A ordem do ranking é frágil, e olhar só o primeiro número é arriscado.
- Um redocking bem-sucedido valida o protocolo para moléculas parecidas naquele sítio. Não valida para outro alvo, como a [Prática B](eif4e.md) mostra.

## Se o servidor estiver fora do ar

Use os arquivos de [resultados/mpro](https://github.com/madsondeluna/curso-docking/tree/main/resultados/mpro): `poses.json` traz energia e RMSD de cada pose, `poses_top3.sdf` traz as três melhores em 3D e `interacoes.json` traz as interações do cristal e das poses. Todos os passos de leitura e comparação continuam possíveis.
