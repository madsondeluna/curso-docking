# Passo a passo: o cap do mRNA na eIF4E

O mesmo protocolo da Prática A, aplicado à eIF4E humana com o cap m7GTP (PDB 1IPC). Aqui ele falha, e a tarefa é descobrir por quê. Os valores vêm das execuções em `resultados/eif4e` e `resultados/eif4e_ac`.

[Voltar à aula](../README.md) · [Passo a passo da Mpro](mpro.md) · [Passo a passo das cavidades](cavidades.md)

## O sistema, e por que ele foi escolhido

Todo mRNA eucariótico começa com um cap: uma guanina metilada, com carga positiva no anel, presa ao resto da fita por três fosfatos. A eIF4E reconhece esse cap e recruta a maquinaria de tradução. É o passo limitante da síntese de proteínas.

**Lado saudável.** Sem esse reconhecimento, a célula não traduz. É o controle de qual mensagem vira proteína.

**Lado patológico.** Células tumorais com eIF4E ativa demais traduzem preferencialmente mensagens ligadas a proliferação e sobrevivência, o que faz da eIF4E um alvo em câncer.

**Por que este caso está na aula.** Ele tem tudo que um exercício bem-comportado não tem: ligante muito carregado, muito flexível e uma interação que a maioria das funções de pontuação trata mal. É assim que a maior parte dos alvos reais se comporta.

**Baixe as entradas:** [dados/eif4e.zip](../dados/eif4e.zip)

## Passo 1: ver o que deveria acontecer

**O que fazer.** Abra [rcsb.org/structure/1IPC](https://www.rcsb.org/structure/1IPC), clique em **3D View** e localize o ligante MGP.

**O que procurar.** A guanina metilada fica **entre dois triptofanos**, Trp56 e Trp102, como o recheio de um sanduíche. Os três fosfatos apontam para fora, em direção a Arg157, Lys162 e Glu103.

**Por que essa geometria é tão específica.** O anel da guanina metilada tem carga positiva, e anéis aromáticos atraem cargas positivas através da nuvem de elétrons π. É a interação cátion-π (Ma e Dougherty, 1997). O sítio da eIF4E é feito para explorá-la: dois anéis, um de cada lado.

**Guarde essa imagem.** Ela é o gabarito contra o qual você vai julgar as poses.

## Passo 2: submeter, com o mesmo protocolo da Prática A

**O que selecionar, campo a campo.**

| Campo | Valor | Por quê |
|-------|-------|---------|
| Ligante | SMILES abaixo, estado de pH 7,4 | mesma regra da Prática A |
| Alvo | `receptor.pdb` do pacote | cadeia A da 1IPC, 1 540 átomos pesados, sem o cap |
| Centro da caixa | 54,59  114,72  −7,80 | centro do cap no cristal |
| Tamanho da caixa | 23 × 22 × 30 Å | extensão do cap mais 8 Å de folga |
| Método | **Docking with AutoDock Vina** | o mesmo da Prática A |
| Exaustividade | 8 | a mesma da Prática A |

```
C[n+]1cn([C@@H]2O[C@H](CO[P@@](=O)([O-])O[P@@](=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]2O)c2nc(N)[nH]c(=O)c21
```

**Repare na carga.** No estado de pH 7,4 o cap tem carga total **−3**, por causa dos três fosfatos. O ensitrelvir da Prática A era neutro. Guarde essa diferença: ela volta no passo 4.

**Por que manter tudo igual à Prática A.** Mudar uma coisa de cada vez é o que permite atribuir a causa. Se o protocolo é o mesmo e o resultado muda, a diferença está no sistema, não no procedimento.

## Passo 3: o resultado, e como lê-lo

| Pose | Energia (kcal/mol) | RMSD (Å) |
|------|--------------------|----------|
| 1 | −6,53 | 6,51 |
| 2 | −6,33 | 6,82 |
| 4 | −6,14 | 4,93 |
| 16 | −4,98 | 4,35 |

**Nenhuma das 20 poses fica abaixo de 3 Å.** A melhor pose em RMSD, a 16, está entre as piores em energia.

**O que olhar além do número.** Abra a pose 1 junto com o cristal. Ela mantém o empilhamento com Trp102, mas **perde o Trp56**: o sanduíche não se forma. O programa achou metade da interação.

**Por que isso é um resultado, e não um fracasso da prática.** Você acabou de ver, com números, o limite do método. Quem só rodou casos fáceis acha que docking resolve; quem viu este caso sabe perguntar o que foi testado.

## Passo 4: levantar e testar hipóteses

**Como proceder.** Liste as causas possíveis, diga o que cada uma prevê e teste o que der para testar.

| Hipótese | Previsão | Teste | Resultado |
|----------|----------|-------|-----------|
| Ligante flexível demais (10 ligações rotáveis) | com menos fosfatos, deve melhorar | repetir com m7GMP | continuou errando: 4,95 Å |
| Função de pontuação fraca em cátion-π e eletrostática | outro método deve melhorar | repetir com Attracting Cavities | melhorou, passo 5 |
| Carga −3 e custo de dessolvatação | ligantes neutros devem ir melhor | comparar com o ensitrelvir | consistente: o neutro acerta |
| Receptor rígido, sítio que fecha sobre o cap | exigiria o receptor em movimento | precisa de dinâmica molecular | fora do alcance desta aula |

**Para testar o m7GMP,** troque a SMILES por esta e repita o passo 2:

```
C[n+]1cn(c2c1C(=O)NC(=N2)N)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)([O-])[O-])O)O
```

**Por que esse teste é informativo mesmo falhando.** Se a flexibilidade fosse a causa principal, tirar dois fosfatos teria resolvido. Não resolveu, o que empurra a explicação para a pontuação e a carga.

## Passo 5: trocar o método, no mesmo servidor

**O que fazer.** Repita a submissão do passo 2 e, no fim, clique em **Docking with Attracting Cavities** em vez do Vina.

**O que muda na tela.** Aparecem parâmetros que só existem nesse método: condições iniciais aleatórias (RIC), exaustividade de amostragem e prioridade de cavidade. A referência usou exaustividade 90, prioridade de cavidade 70 e 2 condições iniciais.

**O que esperar.** Cerca de **15 minutos**, contra 2 do Vina, e 106 poses agrupadas em 50 famílias.

| Medida | Vina | Attracting Cavities |
|--------|------|---------------------|
| Melhor pose do ranking | 6,51 Å | 4,60 Å |
| Menor RMSD entre todas as poses | 4,35 Å | 2,41 Å |
| Poses abaixo de 3 Å | nenhuma | 6% |
| Tempo | 2 min | 15 min |

**O detalhe que importa.** Numa das poses de topo, a base metilada cai praticamente em cima do lugar certo, com RMSD de 0,85 Å considerando só o anel. O que fica fora são os fosfatos.

**A leitura honesta.** O método mais caro encontra a resposta, mas não a coloca em primeiro lugar. Sem o cristal para comparar, você não saberia qual das 106 poses escolher. Trocar de método melhorou a amostragem; não resolveu a pontuação.

## O que concluir

- Um protocolo validado num sistema não está validado em outro. Redocking vale para o alvo em que foi feito.
- Ligantes muito carregados e muito flexíveis são o caso difícil das funções de pontuação rápidas.
- Quando o resultado é ruim, troque uma variável de cada vez e registre o que cada troca previa.
- Sem referência experimental não há como validar uma pose por dentro. É por isso que a próxima etapa natural é a dinâmica molecular, que testa se a pose sobrevive em água.

## Se o servidor estiver fora do ar

Os resultados completos estão em [resultados/eif4e](https://github.com/madsondeluna/curso-docking/tree/main/resultados/eif4e), [resultados/eif4e_ac](https://github.com/madsondeluna/curso-docking/tree/main/resultados/eif4e_ac) e [resultados/eif4e_m7gmp](https://github.com/madsondeluna/curso-docking/tree/main/resultados/eif4e_m7gmp).
