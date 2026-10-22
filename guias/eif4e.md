# Passo a passo: o cap do mRNA na eIF4E

O mesmo protocolo da Prática A, aplicado à eIF4E humana com o cap m7GTP (PDB 1IPC). Aqui ele falha, e o objetivo é descobrir por quê. Os valores de referência vêm das execuções em `resultados/eif4e` e `resultados/eif4e_ac`.

[Voltar à aula](../README.md) · [Passo a passo da Mpro](mpro.md)

## O sistema

A eIF4E reconhece o cap de todo mRNA eucariótico e recruta a maquinaria de tradução. O cap é uma guanina metilada com carga positiva no anel, ligada a três fosfatos.

No cristal, a guanina metilada fica entre Trp56 e Trp102, numa interação cátion-π, e os fosfatos fazem pontes salinas com Arg157, Lys162 e Glu103.

Baixe as entradas: [dados/eif4e.zip](../dados/eif4e.zip).

## 1. Submeter no SwissDock, com Vina

Mesmos passos do guia da Mpro, com estes valores:

**Ligante,** SMILES em pH 7,4, carga total −3:

```
C[n+]1cn([C@@H]2O[C@H](CO[P@@](=O)([O-])O[P@@](=O)([O-])OP(=O)([O-])[O-])[C@@H](O)[C@H]2O)c2nc(N)[nH]c(=O)c21
```

**Alvo:** `receptor.pdb` (cadeia A da 1IPC, 1 540 átomos pesados).

**Caixa:** centro 54,59  114,72  −7,80; tamanho 23 × 22 × 30 Å.

**Método:** AutoDock Vina, exaustividade 8.

## 2. O resultado

| Pose | Energia (kcal/mol) | RMSD (Å) |
|------|--------------------|----------|
| 1 | −6,53 | 6,51 |
| 2 | −6,33 | 6,82 |
| 4 | −6,14 | 4,93 |
| 16 | −4,98 | 4,35 |

Nenhuma das 20 poses fica abaixo de 3 Å. A pose 1 mantém o empilhamento com Trp102, mas perde o Trp56: o sanduíche de triptofanos não se forma.

Compare com a Prática A: a mesma caixa bem posicionada, o mesmo programa, os mesmos parâmetros, e o resultado é outro.

## 3. Investigar a causa

Um resultado ruim não é o fim do exercício. Levante hipóteses e teste as que puder.

| Hipótese | Como testar | O que eu obtive |
|----------|-------------|------------------|
| O ligante é flexível demais (10 ligações rotáveis) | repetir com o m7GMP, que tem só um fosfato | continuou errado: melhor pose a 4,95 Å |
| A função de pontuação não trata bem cátion-π e eletrostática | repetir com o Attracting Cavities, no mesmo servidor | melhorou, ver abaixo |
| A carga −3 cobra muito na dessolvatação | comparar com ligantes neutros, como o ensitrelvir | consistente: o ensitrelvir é neutro e acerta |
| O receptor rígido não reproduz o fechamento do sítio | exigiria dinâmica molecular | não testado aqui |

**Para testar o m7GMP,** use esta SMILES no lugar da anterior:

```
C[n+]1cn(c2c1C(=O)NC(=N2)N)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)([O-])[O-])O)O
```

## 4. Trocar o método no mesmo servidor

Repita a submissão escolhendo **Attracting Cavities** em vez do Vina. Os parâmetros próprios desse método são exaustividade 90, prioridade de cavidade 70 e 2 condições iniciais aleatórias.

**Tempo medido:** 15 minutos, 106 poses em 50 agrupamentos.

| Medida | Vina | Attracting Cavities |
|--------|------|---------------------|
| Melhor pose do ranking | 6,51 Å | 4,60 Å |
| Menor RMSD entre todas as poses | 4,35 Å | 2,41 Å |
| Poses abaixo de 3 Å | nenhuma | 6% |
| Tempo | 2 min | 15 min |

Numa das poses de topo do Attracting Cavities, a base metilada cai praticamente em cima do lugar certo, com RMSD de 0,85 Å considerando só o anel; o que fica fora são os fosfatos.

**A conclusão honesta:** o método mais caro encontra a resposta, mas não a coloca em primeiro lugar. Sem o cristal para comparar, você não saberia qual pose escolher.

## 5. O que concluir

- Um protocolo validado num sistema não está validado em outro.
- Ligantes muito carregados e muito flexíveis são os casos mais difíceis para funções de pontuação rápidas.
- Quando o resultado é ruim, troque uma variável de cada vez: ligante, método, caixa, receptor.
- Sem referência experimental, não há como dizer por dentro se uma pose está certa. É isso que torna o redocking indispensável.

## Se o servidor estiver fora do ar

Os resultados completos estão em [resultados/eif4e](https://github.com/madsondeluna/curso-docking/tree/main/resultados/eif4e), [resultados/eif4e_ac](https://github.com/madsondeluna/curso-docking/tree/main/resultados/eif4e_ac) e [resultados/eif4e_m7gmp](https://github.com/madsondeluna/curso-docking/tree/main/resultados/eif4e_m7gmp).
