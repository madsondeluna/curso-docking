# Passo a passo: encontrar cavidades e decidir quais são sítios

Previsão de cavidades no PrankWeb, que roda o P2Rank, em quatro estruturas escolhidas para responder a uma pergunta cada. Os valores vêm das execuções guardadas em `resultados/cavidades`.

[Voltar à aula](../README.md) · [Passo a passo da Mpro](mpro.md) · [Passo a passo da eIF4E](eif4e.md)

## Por que esta prática vem antes do docking, na vida real

Nas Práticas A e B o sítio já era conhecido: o cristal trazia o ligante no lugar. Fora da aula, você costuma ter só a proteína. Aí a primeira pergunta não é "qual a melhor pose", e sim **"onde eu coloco a caixa"**.

A resposta errada mais comum é supor que o maior buraco da superfície é o sítio. Esta prática mostra, com números, por que isso falha.

## Passo 1: submeter uma estrutura

**O que fazer.** Abra [prankweb.cz](https://prankweb.cz/).

**O que selecionar.**

| Campo | Escolha | Por quê |
|-------|---------|---------|
| Tipo de entrada | **Experimental structure** | vamos usar estruturas resolvidas, não previstas |
| Código | `8HBK` | a Mpro da Prática A |
| Conservação | ligada | resíduos conservados entre espécies costumam ser funcionais, o que ajuda a separar sítio de buraco |

**O que esperar.** Uma fila curta e, em poucos minutos, a página de resultados.

**Por que a conservação ajuda.** Um sítio de ligação sofre pressão evolutiva: mutações ali costumam quebrar a função. Uma reentrância qualquer da superfície não sofre essa pressão. Conservação alta é um indício a mais de que a cavidade faz algo.

## Passo 2: ler a página de resultados

A página tem três áreas: a lista de cavidades, a sequência e o visualizador 3D.

**O que olhar na lista.**

| Coluna | O que significa | Como usar |
|--------|-----------------|-----------|
| Rank | posição pelo escore do P2Rank | ordem de suspeita, não de certeza |
| Score | escore bruto do método | compare entre cavidades da **mesma** proteína |
| Probability | probabilidade calibrada de ser um sítio | leia como confiança, não como verdade |
| Número de resíduos | tamanho da cavidade | cavidades grandes acomodam ligantes maiores |

**O que fazer no visualizador.** Clique na lupa de cada linha para centralizar a cavidade. Os ligantes presentes na estrutura aparecem em bastões: é com eles que você vai comparar.

**O que você deve ver na 8HBK.** Oito cavidades. A primeira, com escore 18,9, coincide com o ensitrelvir. As outras sete ficam de 8 a 26 Å dele.

**A conclusão do primeiro caso.** Uma proteína de 300 resíduos tem oito cavidades e um sítio. **Nem toda cavidade é um sítio.**

## Passo 3: um sítio alostérico aparece?

**O que fazer.** Repita o passo 1 com o código `5MO4`.

**O que essa estrutura tem de especial.** A quinase ABL com **dois fármacos ligados ao mesmo tempo**: o nilotinibe no sítio do ATP e o asciminibe no bolso do miristoil, que é alostérico.

**O que esperar.**

| Cavidade | Escore | Fica a | Que sítio é |
|----------|--------|--------|-------------|
| 1 | 27,9 | 2,1 Å do nilotinibe | ortostérico, o sítio do ATP |
| 2 | 5,2 | 16 Å de qualquer ligante | sem ligante conhecido |
| 3 | 5,0 | 0,9 Å do asciminibe | alostérico, o bolso do miristoil |

**Como interpretar.** O método acha os dois sítios reais, mas põe o alostérico em terceiro lugar, com escore cinco vezes menor que o ortostérico. Sítios alostéricos costumam ser menores e menos óbvios geometricamente.

**Por que isso vale ouro na clínica.** O asciminibe ocupa o lugar onde normalmente entra a cauda miristoilada da própria proteína e trava a quinase numa forma inativa (Wylie *et al.*, 2017). Como não disputa o sítio do ATP, continua funcionando em leucemias com a mutação T315I, que torna a proteína resistente ao imatinibe. Trocar de cavidade foi a saída para a resistência.

**O que você deve tirar daqui.** Se tivesse descartado tudo abaixo do primeiro lugar, teria jogado fora a cavidade que virou fármaco.

## Passo 4: uma cavidade pode não existir ainda

**O que fazer.** Rode `6OIM` e depois `4OBE`. São a mesma proteína, o KRAS: a primeira com o fármaco sotorasibe ligado, a segunda sem fármaco nenhum.

**O que esperar.**

| Estrutura | Cavidade do switch II | Escore |
|-----------|------------------------|--------|
| 6OIM, com o fármaco | aparece em 2º lugar, a 0,9 Å do sotorasibe | 6,7 |
| 4OBE, sem o fármaco | aparece em 3º lugar, rasa | 9,5, contra 21,5 do sítio do GDP |

**A história por trás.** O KRAS foi considerado não medicável por décadas: a superfície é lisa, sem bolsão evidente. Ostrem *et al.* (2013) mostraram que inibidores covalentes da cisteína 12 **abrem** um bolso novo, o switch II, e disso veio o sotorasibe (Canon *et al.*, 2019).

**O que a comparação ensina.** O bolso existe na estrutura sem fármaco, mas raso e mal pontuado. Uma cavidade pequena demais hoje pode ser um sítio amanhã, se algo a abrir. Sítios assim se chamam **crípticos**, e é a dinâmica molecular que costuma encontrá-los.

## Passo 5: julgar uma cavidade

Antes de apontar a caixa do docking para uma cavidade, passe por esta lista.

| Pergunta | Se a resposta for sim |
|----------|------------------------|
| A cavidade tem glicerol, PEG ou sulfato no cristal? | desconfie: pode ser artefato da cristalização |
| Ela fica entre duas cópias da proteína no cristal? | pode não existir na célula |
| Os resíduos são conservados? | indício a favor |
| Ela acomoda uma molécula do tamanho que você quer encaixar? | indício a favor |
| Existe alguma estrutura da mesma proteína com ligante ali? | é a evidência mais forte que existe |
| O escore é alto, mas nada se liga ali em nenhuma estrutura? | escore alto diz que **parece** sítio, não que seja |

**A regra prática.** Escore do P2Rank ordena suspeitas. Só estrutura experimental, ensaio ou mutagênese confirma um sítio.

## O que concluir

- Proteínas têm muito mais cavidades do que sítios: oito contra um, na Mpro.
- Sítios alostéricos são cavidades de verdade, só menores e pior pontuadas; descartá-las por causa do posto é perder alvos.
- Sítios crípticos podem estar rasos ou ausentes na estrutura que você baixou.
- Previsão de cavidade escolhe onde procurar. Ela não substitui evidência de que algo se liga ali.

## Se o servidor estiver fora do ar

Os resultados das quatro estruturas, com todas as cavidades e a distância até cada ligante do cristal, estão em [resultados/cavidades](https://github.com/madsondeluna/curso-docking/tree/main/resultados/cavidades).
