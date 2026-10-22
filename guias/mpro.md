# Passo a passo: o ensitrelvir na Mpro do SARS-CoV-2

Redocking do ensitrelvir na protease principal do SARS-CoV-2 (PDB 8HBK), no SwissDock, com AutoDock Vina. Todos os valores de referência vêm da execução que está em `resultados/mpro`.

[Voltar à aula](../README.md) · [Passo a passo da eIF4E](eif4e.md)

## O que você precisa

Navegador e conexão. Nada é instalado. Nenhum dos serviços pede cadastro.

Baixe o pacote de entradas: [dados/mpro.zip](../dados/mpro.zip). Dentro dele:

| Arquivo | Conteúdo |
|---------|----------|
| `receptor.pdb` | cadeia A da 8HBK, só proteína, sem água nem ligante |
| `ligante_cristal.sdf` | o ensitrelvir na pose do cristal, o gabarito |
| `ligante.smi` | SMILES do ensitrelvir em pH 7,4 |
| `ligante_3d.sdf` e `.mol2` | o ligante em 3D com hidrogênios |
| `caixa.json` | centro e tamanho da caixa de busca |

## 1. Olhar a estrutura antes de tudo

Abra [rcsb.org/structure/8HBK](https://www.rcsb.org/structure/8HBK). Confirme três coisas: método de raios X, resolução 1,8 Å e a presença do ligante 7YY. Clique em **3D View** para abrir o Mol* e localize o ligante dentro da proteína.

Essa checagem é parte do método. Uma estrutura de baixa resolução, ou sem ligante no sítio, muda o que se pode esperar do docking.

## 2. Separar receptor e ligante

O `receptor.pdb` do pacote já está pronto. Ele foi feito mantendo só os átomos `ATOM` da cadeia A, sem água, íons ou o próprio ligante. Para refazer esse passo a partir do PDB bruto, use `tools/prepare_inputs.py`.

Tirar o ligante é o que transforma o exercício num teste: o programa não pode ver a resposta.

## 3. Submeter no SwissDock

Abra [swissdock.ch](https://www.swissdock.ch/) e preencha nesta ordem.

**Ligante.** Cole a SMILES de `ligante.smi`:

```
Cn1cnc(Cn2c(=O)nc(Nc3cc4cn(C)nc4cc3Cl)n(Cc3cc(F)c(F)cc3F)c2=O)n1
```

O servidor mostra os estados de protonação possíveis. Escolha o de pH 7,4. Para o ensitrelvir, a carga total é 0.

**Alvo.** Envie o `receptor.pdb`.

**Método.** Escolha **AutoDock Vina**.

**Caixa.** Use os valores de `caixa.json`:

| Campo | Valor |
|-------|-------|
| Centro X, Y, Z | −21,75  −11,26  2,66 |
| Tamanho X, Y, Z | 28  27  26 |

**Amostragem.** Exaustividade 8.

O servidor estima o tempo antes de enviar e recusa trabalhos longos demais. Se ele recusar, reduza a caixa ou a exaustividade.

**Tempo medido:** 2 minutos, com 20 poses de saída.

## 4. Ler o resultado

Baixe os resultados. O arquivo que interessa é `vina_dock.pdbqt`: cada bloco `MODEL` é uma pose, e a linha `REMARK VINA RESULT` traz a energia em kcal/mol.

Valores que eu obtive:

| Pose | Energia (kcal/mol) | RMSD contra o cristal (Å) |
|------|--------------------|----------------------------|
| 1 | −9,13 | 1,51 |
| 2 | −8,71 | 8,32 |
| 3 | −8,64 | 6,73 |
| 4 | −8,43 | 7,27 |

Suas energias podem variar um pouco: o Vina usa uma semente aleatória. Espere a pose 1 entre −9,5 e −8,5 kcal/mol, e abaixo de 2 Å.

## 5. Calcular o RMSD

O RMSD compara a pose com o ligante do cristal, átomo pesado a átomo pesado, sem realinhar as estruturas: receptor e cristal já estão no mesmo referencial.

Sem instalar nada, dá para fazer a comparação visual: abra `ligante_cristal.sdf` e a pose 1 no mesmo visualizador e veja se os anéis caem no mesmo lugar. Para o número, `tools/analyze_poses.py` faz o cálculo considerando a simetria do ligante.

Critério da área: abaixo de 2 Å a pose está certa, entre 2 e 3 Å está perto, acima de 3 Å está errada.

## 6. Ver as interações no PLIP

Monte um arquivo PDB com o receptor e a pose juntos e envie para o [PLIP](https://plip-tool.biotec.tu-dresden.de/). Ou, mais rápido, envie o código 8HBK e veja as interações do cristal, para comparar com o que a sua pose faz.

Interações do cristal, pelo PLIP:

| Tipo | Resíduos |
|------|----------|
| Ligação de hidrogênio | Thr26, Gly143, Ser144, Cys145, Glu166 |
| Hidrofóbica | Thr25, Met165 |
| Empilhamento π | His41 |

A pose 1 que eu obtive repete 75% dessas interações: perde o contato com Thr25 e Thr26, na borda do sítio, e mantém todo o resto.

## 7. O que concluir

- A pose de melhor energia é a pose certa. O protocolo funciona neste sistema.
- A diferença de energia entre a pose 1 e a pose 2 é de 0,4 kcal/mol, e a pose 2 está 8 Å fora. A ordem do ranking é frágil.
- Um redocking que dá certo valida o protocolo para moléculas parecidas no mesmo sítio. Não valida para outro alvo.

## Se o servidor estiver fora do ar

Use os arquivos de [resultados/mpro](https://github.com/madsondeluna/curso-docking/tree/main/resultados/mpro): `poses.json` traz energia e RMSD de cada pose, `poses_top3.sdf` traz as três melhores em 3D e `interacoes.json` traz as interações do cristal e das poses.
