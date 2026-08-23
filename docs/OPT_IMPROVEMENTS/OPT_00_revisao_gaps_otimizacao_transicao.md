# Otimização Aerodinâmica com Transição: Estado da Arte e Gaps para Contribuição de Doutoramento

Objetivo final da tese: um método de otimização adjunta com transição que produza resultados robustos em configurações completas de avião e que seja capaz de informar decisões de engenharia em processos de design de alta cadência.

---

## 1. Estado da arte por eixo

### 1.1 Modelos de transição compatíveis com RANS
- Langtry e Menter (2009, AIAA J): modelo γ-Reθt local baseado em correlações, tornou a transição viável em códigos paralelos não estruturados. Tese de Langtry (2006) documenta as correlações.
- Menter et al. (2015): modelo γ de uma equação, Reθt passa a algébrico e local. Variante DLR-γ de Grabe e Krumbein com extensões de crossflow.
- Coder e Maughmer (2014, AIAA J) e Coder (2019): Amplification Factor Transport, AFT2019b inclui crossflow. Base física em e^N em vez de correlações de intermitência.
- Piotrowski e Zingg (2021 a 2023, AIAA J e teses): variante suave SA-sLM2015, restrição de passo de tempo nos termos fonte para solver totalmente implícito, correções de compressibilidade para TS e crossflow (sLM2015cc). É a tua base de implementação.
- Righi, Halila e Azevedo (2026): suavização do termo de crossflow, confirma que a suavização de termos fonte continua a ser tema aberto.

### 1.2 Otimização adjunta com transição
Linha temporal essencial:
- Driver e Zingg (2007): primeiro acoplamento e^N com adjunto Newton-Krylov em perfis.
- Lee e Jameson (2009): módulo de transição com bases e^N duplas (TS e CF) acoplado a adjunto contínuo, perfis e asas, mas com acoplamento fraco e turbulência algébrica.
- Khayatzadeh e Nadarajah (2011 a 2014): adjunto discreto do γ-Reθt para perfis NLF.
- Rashad e Zingg (2016): e^N fortemente acoplado no RANS e no gradiente, otimização multiponto de perfis NLF.
- Shi et al. (2020, AIAA J): adjunto discreto de SA-AFT em perfis, no framework MACH-Aero.
- Halila, Fidkowski e Martins (2019 a 2021): adjunto RANS transicional, análise de sensibilidades e tese sobre transição em otimização de alta fidelidade.
- Djeddi e Ekici (2022, 2024, AST): AFT-S suave, e primeiro trabalho a usar diretamente o atraso da posição de transição como objetivo em 2D, com deteção diferenciável da posição via sigmoide do índice de turbulência. Precedente direto da tua ideia do funcional da frente, mas apenas em perfis.
- Paredes et al. (2024, NASA/AIAA): adjunto do problema de valores próprios da LST (dual N-factor TS e CF) acoplado ao FUN3D, perfis com e sem flecha.
- Fizaa e Zingg (2024): otimização de asas em flecha (infinitas e finitas) com transição livre, com restrições geométricas para suprimir attachment-line e Görtler.
- Shi et al. (2018, AIAA J): alternativa surrogate (kriging) para asa NLF transónica com e^N duplo, 42 variáveis, sem adjunto.

Leitura do estado da arte: o adjunto com transição está maduro em 2D e a chegar a asas simples. Ninguém demonstrou otimização transicional robusta de configuração completa (asa, fuselagem, winglet, propulsão).

### 1.3 Robustez e incerteza
- Zhao et al. (2023, Chinese J of Aeronautics): otimização robusta de asa laminar com adjunto e polynomial chaos com gradiente, incerteza em condições de voo, base MACH-Aero com e^N.
- Rashad e Zingg (2016): multiponto como proxy de robustez.
- Gap conhecido na comunidade: sensibilidade dos ótimos NLF a Tu_∞ e N_crit raramente quantificada com modelos LCTM, e nunca com formulação de margens.

### 1.4 Malha e adaptação
- Requisitos de malha dos modelos LCTM (documentados desde Langtry 2006): y+ da ordem de 1 e resolução em corda muito fina na zona da frente. Nenhum estudo quantifica o erro do gradiente adjunto em função da resolução na frente.
- Adaptação orientada a outputs (Fidkowski e Darmofal, 2011, revisão) existe para RANS, e o grupo de Fidkowski aplicou adaptação a escoamentos transicionais em contexto de análise. Não existe adaptação dentro do loop de otimização com tratamento da consistência do gradiente quando a malha muda.

### 1.5 Modos não simulados e a sua mitigação
- Attachment-line contamination: critério de Poll (1979) com R̄ crítico. Trabalhos experimentais e DNS (Spalart 1989) estabelecem os limiares. Na otimização, o estado da arte é tratá-lo com restrições geométricas (Fizaa e Zingg 2024) ou dispositivos (fences, Gaster bumps). Não existe critério de contaminação diferenciável avaliado a partir da solução CFD e imposto como restrição adjunta.
- Görtler, bypass, rugosidade e steps: tratados por margens empíricas de ΔN (por exemplo Crouch para steps) fora do loop de otimização. Nunca dentro.
- CATNLF (Campbell e Lynde, NASA, e trabalhos Bayesianos recentes em 2026): filosofia de pressão alvo para atenuar crossflow a flecha alta, via design inverso, não via adjunto com modelo LCTM.

### 1.6 Efeito de hélice
- Avaliação de modelos de transição em slipstream de propulsão distribuída (AIAA 2025): DLR-γ concorda com e^N, γ-Reθt transiciona cedo demais, e em regime instacionário os vórtices de ponta de pá deslocam ciclicamente a frente.
- Estudos experimentais clássicos: tractor antecipa a transição, pusher pode estendê-la.
- Não existe nenhuma otimização adjunta com transição sob slipstream. Nem sequer com actuator disk estacionário.

### 1.7 IA em transição
- Substituição da LST por redes: CNN treinada em perfis de camada limite para prever taxas de amplificação (Zafar, Xiao et al. 2020), modelos de previsão de posição de transição por NN em várias configurações, field inversion and machine learning (Duraisamy) para corrigir modelos.
- Nada disto foi construído de forma analiticamente diferenciável para entrar num adjunto discreto de otimização. E nenhum surrogate de frente é usado para guiar malha ou otimizador.

---

## 2. Gaps e contribuições propostas

Cada gap tem novidade defensável para PhD. Ordenados pelos teus eixos.

### Solver (análise direta mais rápida e robusta)
**G1. Condicionamento do sistema transicional.** Reescalonamento dos blocos de estado, tratamento implícito exato do fonte de relaxação do Reθt, precondicionador field-split PETSc por bloco físico, e deflação ou reciclagem de modos lentos no Krylov ao longo da otimização. Não há tratamento sistemático publicado para Newton-Krylov com LCTM em geometrias complexas. Entrega: convergir onde hoje falha, com o mesmo input.

### Equações
**G2. Modelo de 2 equações no total (SA + γ, Reθt algébrico), suave e compressível.** O análogo SST-γ existe, o análogo SA suave com correções de compressibilidade e adjunto 3D não. Menos um campo transportado é menos rigidez, menos memória de adjunto e menos a equação que te trava.
**G3. Regularização da espessura da frente com continuação.** Parametrizar a largura da zona de transição e fazer continuação do difuso para o afiado. Não existe na literatura de transição, tem análogo forte em topologia e phase-field, e ataca diretamente a exigência de malha.

### Adjunto
**G4. Adjunto com operador suavizado.** Primal afiado, adjunto do Jacobiano regularizado como direção de descida, com região de confiança para gradientes inexatos e fase final com adjunto exato. Com o teu adjunto verificado consegues medir o ângulo entre gradiente aproximado e exato, coisa que ninguém publicou.
**G5. Erro do gradiente versus resolução na frente.** Estudo quantitativo com critério de resolução derivado (Δx sobre comprimento de transição). Fundação citável para tudo o resto e não existe.

### Malha
**G6. Adaptação anisotrópica da banda da frente dentro do loop de otimização.** Indicador construído com o adjunto de γ, esquema de adaptação congelada por ciclo com verificação de consistência do gradiente. Variante com IA: surrogate que prevê a frente a partir do Cp e da geometria para pré-adaptar a malha antes de cada solve.

### Processo de otimização
**G7. Funcional diferenciável da posição da frente em 3D.** Djeddi e Ekici fizeram o precedente 2D. Aberto: formulação em asa (centróide do gradiente de γ por faixa de envergadura), restrição de regularidade da frente em envergadura, e estudo de condicionamento de Cd puro contra objetivo composto, que explica a tua estagnação.
**G8. Robustez com margens.** Otimização com momentos por adjunto sobre Tu_∞ e N_crit com modelo LCTM (o precedente de 2023 usa e^N), e formulação de margem de transição em vez de valor esperado.

### Modos não simulados, restrições inteligentes
**G9. Restrições diferenciáveis de margem para mecanismos não resolvidos.** Critério de Poll avaliado na solução CFD ao longo da linha de atrelamento (raiz e junção asa-winglet), margem de Görtler, penalização de sensibilidade a rugosidade e steps via ΔN, tudo agregado com KS e imposto no SNOPT. Transforma "o modelo não simula isto" em "o design garante margem contra isto". É a peça que torna os resultados críveis para engenharia e não existe em nenhum framework adjunto.

### Efeito de hélice
**G10. Otimização transicional sob slipstream.** Primeiro nível: actuator disk estacionário acoplado ao teu modelo, com Tu_∞ local elevado no tubo de corrente como campo de contaminação, e otimizar a asa com e sem hélice para quantificar a perda de laminaridade recuperável por design. Segundo nível: modelo estacionário de contaminação cíclica (ponderação do estado laminar e turbulento pelo duty cycle das esteiras das pás, na linha dos modelos de wake-induced transition de turbomáquinas). Nenhuma das duas coisas foi feita com adjunto.

### Aplicação integradora
**G11. Configuração completa.** Asa, fuselagem com contaminação na junção, winglet com zona não linear de crossflow, hélice. É a demonstração que junta tudo e que ninguém tem.

---

## 3. Roadmap de artigos

1. **P1, fundação numérica:** G1 e G5. Condicionamento e critério de resolução do gradiente. Journal alvo tipo AIAA J ou JCP conforme o peso teórico.
2. **P2, método:** G3 e G4. Continuação da frente e adjunto suavizado, com a métrica do ângulo do gradiente. Contribuição central da tese.
3. **P3, modelo:** G2. SA + γ com Reθt algébrico, suave, compressível, com adjunto 3D validado nos teus casos (placa, NLF(1)-0416, S809, NLF(2)-0415, Sickle Wing).
4. **P4, otimização:** G7 e G9. Funcional da frente em 3D mais restrições de margem (attachment line na junção). Primeira otimização com garantias contra modos não simulados.
5. **P5, aplicação:** G8, G10 e G11. Asa completa com winglet e slipstream, com robustez. O paper que fecha a frase da tese.

## 4. Riscos de colisão a vigiar
- Djeddi e Ekici podem estender o transition delay a 3D. Diferencia-te pela regularidade em envergadura, condicionamento e margens.
- Grupo de Zingg (Fizaa) está ativo em asas em flecha com restrições para attachment line. Diferencia-te pelo critério flow-based diferenciável em vez de proxy geométrico.
- Grupos chineses (Shi, Zhao) dominam robustez com e^N. Diferencia-te pelo LCTM adjunto e pelas margens.

---

# ADENDA (mantém tudo acima, isto acrescenta)

## 5. SA-BCM como peça estratégica

O que existe na literatura:
- Cakmakcioglu, Bas et al. (2017) e a revisão BCM (Mura e Cakmakcioglu, AIAA 2020-2714): modelo de intermitência algébrica que apenas amortece a produção do SA, zero equações extra, implementado no SU2.
- Existe já um framework de otimização adjunta NLF com SA-BC no SU2 (AIAA Journal, 2021, doi 10.2514/1.J059923): adjunto discreto por AD, validado contra diferenças finitas, aplicado ao NLF(1)-0416 e a uma asa de baixo alongamento. Ou seja, o adjunto BCM simples já foi feito. A novidade não está em repetir isso, está em usar o BCM como ferramenta dentro de um framework de fidelidade superior.

Quatro usos publicáveis do teu BCM já implementado:

**5.1 Multi-fidelidade com gestão de modelos.** Trust region model management (Alexandrov et al., 2001) com BCM como fidelidade baixa e sLM2015 como alta. O TRMM clássico usa Euler vs RANS ou malha grossa vs fina. Nunca foi feito com dois modelos de transição de custo e rigidez diferentes. O BCM explora, o sLM2015 corrige e certifica.

**5.2 Adjunto de modelo trocado.** Primal com sLM2015, adjunto montado com o Jacobiano do BCM. É o análogo transicional do frozen turbulence adjoint, mas com um modelo de transição consistente em vez de turbulência congelada. Como tens o adjunto sLM2015 exato e verificado, consegues medir o ângulo entre a direção BCM e o gradiente verdadeiro ao longo de uma otimização real. Ninguém publicou essa quantificação para transição. Liga diretamente ao G4 e é mais barato que o operador suavizado, porque o BCM já existe.

**5.3 Continuação entre modelos.** Fase 1 da otimização com BCM (suave, sem Reθt, sem os teus problemas de condicionamento), fase 2 arranca do design da fase 1 com sLM2015. Homotopia em fidelidade de modelo em vez de homotopia em parâmetro. Custo quase nulo de implementação e resultado imediato para o teu caso de estagnação.

**5.4 Estudo de condicionamento head-to-head.** Espectro e número de condição do Jacobiano BCM vs sLM2015 na mesma malha e no mesmo caso. Isola quantitativamente de onde vem a rigidez (Reθt transportado, correlações, acoplamento). É a evidência empírica que sustenta o argumento do modelo de 2 equações (G2) e serve de secção forte no P1.

Limitações a declarar sempre: o BCM não tem crossflow, recebe Tu como input direto e tem menos física. É motor de exploração, não de certificação.

## 6. Ideias fora da caixa, importadas de outros campos

**6.1 Phase-field e Γ-convergência (matemática da fratura e topologia).** Tratar γ como campo de fase tipo Allen-Cahn com largura de interface ε controlada. A teoria de Γ-convergência (Bourdin e Chambolle para topologia phase-field) dá o enquadramento para afirmar que os ótimos do problema regularizado convergem para os do problema afiado quando ε tende ao valor físico. É a espinha teórica que transforma o G3 de truque numérico em contribuição matemática. Encaixa no teu interesse já declarado em phase-field e SAV.

**6.2 Robustez erode/dilate (otimização topológica, Wang, Lazarov e Sigmund, 2011, SMO).** Em topologia, robustez a erros de fabrico obtém-se otimizando em simultâneo o design erodido, o intermédio e o dilatado. Tradução direta: otimizar em simultâneo as frentes com N_crit menos ΔN, N_crit e N_crit mais ΔN. Três primals por iteração, formulação minimax ou soma pesada, e o resultado é um design cuja laminaridade não colapsa com a receptividade real do voo. Barato, elegante, e resolve o teu problema de sensibilidade ao input com uma ideia que ninguém trouxe para transição.

**6.3 Transition fitting (analogia com o revival do shock-fitting, Paciorri e Bonfiglioli, 2009 em diante).** Em vez de capturar a frente na malha, tratá-la como incógnita explícita do problema: uma interface com condição de equilíbrio (por exemplo N igual a N_crit ou o critério de onset do modelo) imposta como equação adicional, com a malha alinhada e deformada com a frente pelo teu pyGeo/IDWarp. O Newton resolve a posição da frente junto com o escoamento. Elimina por construção a exigência de malha fina, porque a descontinuidade vive entre células e não dentro delas. Alto risco, prémio máximo: seria uma mudança de paradigma na forma de simular e otimizar transição.

**6.4 Fronteira livre e cálculo de formas (problemas de Stefan).** A frente de transição como fronteira livre com derivada de Hadamard bem definida da sua posição em relação à forma. A matemática de adjuntos de problemas de fronteira livre existe (otimização de mudança de fase, obstáculos). Dá um formalismo alternativo e rigoroso ao funcional da frente do G7, com teoremas em vez de heurísticas.

**6.5 Harmonic balance com transição para o efeito de hélice.** O adjunto harmonic balance é maduro em turbomáquinas (Huang e Ekici, 2014, AST) e existe no SU2 para escoamentos quase periódicos (Rubino et al., 2018, JCP). Modelos de transição dentro de HB não existem, e adjunto de HB transicional muito menos. A contaminação cíclica das esteiras das pás é exatamente um fenómeno periódico de frequência conhecida, o caso ideal para HB. Otimizar a frente de transição periódica sob slipstream seria a primeira otimização transicional instacionária da literatura. Flagship natural para o G10 nível 2.

**6.6 Neural operators diferenciáveis dentro do adjunto.** O precedente CNN-LST (Zafar e Xiao, 2020) prevê taxas de amplificação a partir de perfis de camada limite, mas sem adjunto. Substituir o bloco de correlações (ou a própria LST) por um operador neural pequeno com derivadas exatas por autodiff, acoplado ao teu adjunto discreto, dá gradientes suaves por construção e remove os smoothMinMax. Ligação natural a field inversion and machine learning (Duraisamy et al., 2019, Annu. Rev. Fluid Mech.) para calibrar com dados experimentais ou de voo (termografia IV), na direção de um gémeo digital laminar.

**6.7 One-shot / SAND (Griewank, Gauger).** Resolver primal, adjunto e atualização de design em simultâneo (piggyback), em vez do ciclo sequencial converge-tudo-depois-otimiza. O ganho é maior precisamente quando o primal é caro e rígido, que é o teu caso. Nunca aplicado a escoamentos transicionais. Nota: Djeddi e Ekici (2021, AIAA J) têm one-shot para instacionário, mais uma razão para vigiar esse grupo.

**6.8 Chance constraints (programação estocástica).** Formalizar as margens do G9 como restrições probabilísticas: probabilidade da frente recuar além de x inferior a ε, avaliada por momentos calculados via adjunto. Transforma "margem de engenharia" em objeto matemático otimizável e liga o G8 ao G9 numa formulação única.

**6.9 Semismooth Newton (otimização não suave).** Em vez de suavizar todos os min/max das correlações (caminho Piotrowski), tratá-los como condições de complementaridade e resolver com métodos semismooth (Hintermüller, Ito, Kunisch). Mantém o modelo original exato e ataca a não-suavidade no solver em vez de a esconder no modelo. Diferenciador teórico direto face à linha UTIAS, e já está na tua lista de interesses.

**6.10 SAV e integradores stiff (gradient flows).** Scalar auxiliary variable e esquemas Rosenbrock para os termos fonte rígidos de γ, importados da comunidade de phase-field (Shen et al., 2018). Estabilidade incondicional demonstrável para a parte reativa do sistema. Contribuição de solver com teoria por trás, complementar ao G1.

## 7. O caminho mais substancial (recomendação focada no teu problema)

Se o critério é resolver o teu problema concreto (frente fina, malha proibitiva, estagnação, avião completo) com máxima publicabilidade, a espinha é esta:

1. **Fundação (P1):** G5 (erro de gradiente vs resolução na frente) mais 5.4 (condicionamento BCM vs sLM2015). Estabelece os números que justificam tudo o resto e já usa só infraestrutura que tens.
2. **Núcleo metodológico (P2, a contribuição central):** regularização phase-field da frente com continuação (G3 com a teoria de 6.1) combinada com o adjunto de operador trocado ou suavizado (G4 com 5.2), tudo medido com a métrica do ângulo do gradiente que só tu consegues calcular. Resolve diretamente a frente fina e a estagnação.
3. **Credibilidade de engenharia (P4):** margens diferenciáveis (G9) com robustez erode/dilate (6.2) e, se quiseres formalizar, chance constraints (6.8). É o que sustenta a frase "capaz de informar decisões de engenharia".
4. **Flagship de aplicação (P5):** junção e winglet com critério de Poll, e hélice. Nível 1 com actuator disk estacionário (G10), nível 2 com HB transicional (6.5) se o tempo der. O nível 2 sozinho vale um artigo de grande visibilidade.

O transition fitting (6.3) fica como aposta paralela de alto risco: se funcionar num caso 2D, é um artigo de destaque por si só e pode reordenar a tese inteira.

## 8. Referências adicionais desta adenda

- Cakmakcioglu, S. C., Bas, O. et al., correlações BCM, AIAA 2017 e AIAA Paper 2020-2706 / 2020-2714 (modelo revisto)
- Discrete Adjoint-Based ASO Framework for Natural Laminar Flows (SA-BC no SU2), AIAA Journal, 2021, doi 10.2514/1.J059923
- Alexandrov, N. et al., trust region model management, 2001
- Wang, F., Lazarov, B., Sigmund, O., robust topology optimization (erode/dilate), Struct. Multidisc. Optim., 2011
- Bourdin, B., Chambolle, A., phase-field em otimização de design
- Paciorri, R., Bonfiglioli, A., unstructured shock-fitting, JCP, 2009 em diante
- Huang, H., Ekici, K., discrete adjoint harmonic balance, Aerospace Science and Technology, 2014
- Rubino, A. et al., adjoint harmonic balance para escoamentos quase periódicos, J. Comput. Phys., 2018
- Djeddi, R., Ekici, K., one-shot para otimização instacionária, AIAA Journal, 2021
- Zafar, M. I., Xiao, H. et al., CNN baseada em LST para transição, 2020
- Duraisamy, K., Iaccarino, G., Xiao, H., turbulence modeling in the age of data, Annu. Rev. Fluid Mech., 2019
- Hintermüller, M., Ito, K., Kunisch, K., semismooth Newton, SIAM J. Optim., 2002
- Shen, J. et al., scalar auxiliary variable para gradient flows, J. Comput. Phys., 2018
- Griewank, A., Gauger, N., one-shot / piggyback optimization

---

# PARTE II: A DENSIDADE DE MALHA (fio prioritário, sem apagar os restantes)

Nota de enquadramento: esta parte aprofunda UM dos fios da tese. Todos os outros (G1 a G11, adenda das secções 5 a 7) mantêm-se em análise. A frase de fecho da tese continua a precisar deles todos.


Foco: resolver a exigência de malha dos modelos de transição em otimização adjunta, em malhas estruturadas multibloco (ADflow), com o alvo final em avião completo.

---

## 1. O problema em números

- A frente de transição tem comprimento l_tr que encolhe com o Reynolds. A exigência de resolução em corda cresce à medida que l_tr encolhe. Isto está documentado pelo próprio Piotrowski: a capacidade do framework de explorar o espaço de design NLF é sensível à resolução em corda nas zonas de transição, e simulações com transição livre precisam de mais malha do que fully turbulent para a mesma precisão.
- O teu upwind de 1ª ordem em γ e Reθt espalha a frente, o que mascara o problema em malha grossa e o devolve em malha fina.
- A tua observação central, que a derivada da frente exige malha ainda mais fina que a análise e que a frente "não se mexe" em malha grossa, não está quantificada em lado nenhum da literatura. Isso é o teu ponto de entrada.

## 2. O que já existe (para não reinventar a roda)

### 2.1 Requisitos de malha para modelos de transição (análise, não otimização)
- Langtry (tese 2006) e NASA TMR: y+ máximo da ordem de 1, recomendações qualitativas de espaçamento em corda e razão de expansão.
- Malan et al. (AIAA 2009-1142, calibração do γ-Reθ): sensibilidade forte ao espaçamento em corda perto do fim da transição, verificada por duplicação de malha.
- Piotrowski e Zingg: sensibilidade da exploração do espaço de design à resolução em corda na frente, requisito cresce com Re. É o mais próximo que existe do teu tema, mas fica pela constatação. Nenhum critério quantitativo, e nada sobre o erro do gradiente.

### 2.2 Prática fully turbulent em otimização (a tua baseline de comparação)
- Prática padrão MACH-Aero e ADODG: famílias de malha L0/L1/L2, otimizar em L1/L2, verificar em L0, Richardson para o contínuo. Estudos de referência (por exemplo Lyu, Kenway e Martins no CRM) documentam a dependência do ótimo em relação ao nível de malha e o fenómeno do otimizador a explorar o erro de discretização (o design "ganha" drag numérico que desaparece na malha fina).
- Warm start multinível (otimizar grosso, refinar, continuar) é usado como aceleração, mas sem controlo de erro formal.
- Conclusão útil: mesmo em fully turbulent o problema da malha em otimização só é tratado por protocolo, não por método. Com transição o protocolo rebenta, porque o requisito local na frente é muito mais severo e a frente muda de sítio durante a otimização.

### 2.3 Adaptação orientada a outputs dentro da otimização (o estado da arte real do teu tema)
A linha Chen e Fidkowski (Michigan) é a referência direta e tens de a citar e distinguir:
- Chen e Fidkowski (AIAA 2017-3102): otimização de perfil com malhas adaptativas por estimativas de erro adjuntas, incluindo erro das restrições.
- Chen e Fidkowski (JCP 2019, 387:163-185): controlo do erro de discretização em otimização com restrições.
- Chen e Fidkowski (AST 2020, doi 10.1016/j.ast.2020.106004): otimização multiponto de fidelidade variável com malhas adaptadas por outputs, começa em malha grossa e refina à medida que o design converge.
- Chen e Fidkowski (Computers and Fluids 2021): estimativa de erro e adaptação por CNN.
Limitações de toda esta linha, e é aqui que entras: discretização DG de alta ordem, malhas não estruturadas, 2D ou quase, fully turbulent ou laminar. Nada de transição, nada de estruturado multibloco, nada de avião completo, e a consistência do gradiente entre ciclos de adaptação é tratada de forma pragmática.
- Antecedentes citados por eles: Nemec e Aftosmis (cart3D, otimização progressiva com refinamento), Dalle e Fidkowski (multi-fidelidade), Li e Hartmann (adaptação multi-alvo).

### 2.4 Adaptação métrica anisotrópica (análise de alta fidelidade)
- Alauzet e Loseille (revisões 2016), Park et al. (refine, NASA), ciclo de ponto fixo malha-solução. Maturidade alta para análise em não estruturado, incluindo casos AIAA workshops. Não vive dentro de loops de otimização com adjunto de design, e não é estruturado.
- Adaptação para escoamentos transicionais em análise: trabalho inicial existe na órbita de Fidkowski e Halila (nota a ligação: Halila, Chen e Shi partilham publicação em AST 2019). É a confluência a vigiar.

### 2.5 r-adaptação e redistribuição (o teu cavalo de batalha natural)
- Huang e Russell (livro Adaptive Moving Mesh Methods), Budd, Huang e Russell (Acta Numerica 2009): teoria de equidistribuição e MMPDEs.
- Budd e Williams (SISC 2009), Browne, Budd et al. (JCP 2014, r-adaptação 3D rápida): redistribuição por Monge-Ampère parabólica, malhas sem emaranhamento, conectividade preservada, usada operacionalmente em previsão meteorológica.
- Propriedade chave para ti: r-adaptação NÃO muda topologia nem conectividade. Numa malha estruturada multibloco isso significa que os blocos, as faces e todo o teu adjunto ficam intactos. É exatamente o tipo de adaptação compatível com o ADflow, e nunca foi aplicada nem a frentes de transição nem dentro de otimização aerodinâmica.

### 2.6 ML para adaptação
- Chen e Fidkowski (2021) com CNN para o estimador de erro, e trabalhos de RL para geração de malha. Embrionário, espaço aberto.

## 3. O gap exato

Cada par já existe em parte. A interseção tripla está vazia:

{transição} ∩ {malha adaptativa} : só constatações de sensibilidade, nenhum método
{malha adaptativa} ∩ {otimização adjunta} : Chen e Fidkowski, mas DG não estruturado, sem transição
{transição} ∩ {otimização adjunta} : a tua área, mas toda em malha fixa
{qualquer coisa} ∩ {estruturado multibloco r-adaptativo} ∩ {avião completo} : vazio

E a pergunta que ninguém fez: qual é o erro do GRADIENTE, não do output, em função da resolução na frente, e o que isso faz ao otimizador. Tu tens a ferramenta única para responder (adjunto exato verificado a 3 vias).

## 4. O que podes fazer, por fases

### Fase 0, o estudo que funda tudo (publicável já)
Família de malhas onde SÓ o espaçamento em corda na zona da frente varia, tudo o resto fixo. Para cada nível: erro de Cd, erro da posição da frente, e erro do gradiente adjunto (norma e ângulo contra a malha mais fina). Três curvas em função de Δx/l_tr. Extras que ninguém tem:
- deteta o "front pinning": a histerese da frente quando o otimizador tenta movê-la em malha grossa (varre um DV, mede a resposta da frente, mostra o degrau)
- deriva o critério: N células por comprimento de transição para gradiente fiável, com N medido e não chutado
- repete em fully turbulent para mostrar o contraste com a prática L0/L1/L2
Resultado: o primeiro critério de resolução para gradientes em transição. É a tua secção de motivação para todos os papers seguintes.

### Fase 1, r-adaptação diferenciável na malha estruturada (a contribuição central)
Redistribuição de pontos ao longo das linhas de malha em corda (equidistribuição 1D por linha, trivial em estruturado) ou Monge-Ampère parabólica para o caso geral, com função monitor construída de ∇γ ou, melhor, do indicador adjunto (sensibilidade de γ pesada pelo adjunto do objetivo).
O argumento matador, específico do teu framework: no MACH-Aero as coordenadas da malha já são input diferenciado da cadeia adjunta (dR/dX existe e é usado pelo warping). Se a redistribuição for um mapa suave e diferenciável do estado e do design, ela entra na cadeia total e o gradiente continua exato COM a adaptação ligada. Ou seja, adaptação de malha totalmente diferenciável dentro do adjunto. A linha Chen e Fidkowski não tem isto, eles re-adaptam e aceitam a inconsistência entre ciclos. Tu podes ter o primeiro esquema de adaptação com consistência de gradiente provada.
Protocolo alternativo mais simples para começar: adaptação congelada por ciclo major do SNOPT com verificação de consistência e restart de região de confiança.

### Fase 2, casar com a regularização da frente
Liga a largura da frente ao espaçamento local: escolhe o parâmetro de suavização para a frente ocupar sempre N células (o N da Fase 0), e vai afiando à medida que a r-adaptação concentra pontos. Regularização ciente da resolução. É a síntese das duas linhas (continuação e adaptação) num único método com garantia: frente sempre resolvida, em qualquer malha, com gradiente sempre informativo.

### Fase 3, avião completo
- Em estruturado multibloco puro, a r-adaptação continua a valer porque não mexe na topologia dos blocos.
- O ADflow tem overset: um patch overset fino que segue a frente (translada e deforma com ela) é uma segunda via de refinamento local sem tocar na malha de fundo. Refinamento local por componente (asa, winglet, junção) sem remeshing global.
- Fecha com o caso demonstrador: asa com winglet, frente resolvida por adaptação, gradiente consistente, custo comparado contra a malha uniformemente fina equivalente.

## 5. Como se posiciona a novidade (frases para a tese)

- Primeiro critério quantitativo de resolução para GRADIENTES em escoamentos transicionais, incluindo o fenómeno de front pinning na otimização.
- Primeira adaptação de malha para transição dentro de um loop de otimização adjunta.
- Primeira r-adaptação em malha estruturada multibloco com consistência de gradiente por diferenciação do mapa de redistribuição.
- Demonstração em configuração com winglet e junção, no caminho do avião completo, com redução de X vezes nos pontos de malha para o mesmo erro de gradiente.

## 6. Riscos e vigilância

- Michigan (Fidkowski, com a ponte Halila) é quem mais facilmente estende adaptação por outputs a transição. A tua defesa é o estruturado multibloco, a diferenciabilidade do mapa de adaptação e o foco no gradiente.
- A equidistribuição 1D por linha pode degradar ortogonalidade perto da parede em geometrias curvas. Mitiga com redistribuição só acima da subcamada ou com o termo de regularização do MMPDE.
- O critério da Fase 0 depende do esquema (o teu upwind de 1ª ordem espalha a frente). Reporta o critério para 1ª ordem e, se fizeres a convecção de 2ª ordem de γ, mostra como o critério muda. Isso por si só é mais uma secção.

## 7. Referências desta nota

- Langtry, R., tese de doutoramento, Univ. Stuttgart, 2006
- Malan, P., Suluksna, K., Juntasaro, E., calibração do γ-Reθ, AIAA 2009-1142
- NASA Turbulence Modeling Resource, páginas dos modelos SA-noft2-Gamma-Retheta
- Piotrowski, M., Zingg, D., artigos AIAA J 2021 a 2023 e tese UTIAS (sensibilidade à resolução em corda na frente)
- Lyu, Z., Kenway, G., Martins, J., benchmark ADODG do CRM, AIAA J 2015 (dependência de malha do ótimo)
- Chen, G., Fidkowski, K., AIAA 2017-3102, JCP 2019 (387:163-185), AST 2020 (106004), Computers and Fluids 2021
- Nemec, M., Aftosmis, M., otimização progressiva com refinamento adaptativo (Cart3D)
- Fidkowski, K., Darmofal, D., revisão de estimativa de erro por outputs, AIAA J, 2011
- Alauzet, F., Loseille, A., revisão de adaptação anisotrópica, 2016. Park, M. et al., refine (NASA)
- Huang, W., Russell, R., Adaptive Moving Mesh Methods, Springer, 2011. Budd, Huang, Russell, Acta Numerica, 2009
- Budd, C., Williams, J., SIAM J. Sci. Comput., 2009 (Monge-Ampère parabólica). Browne, P., Budd, C. et al., JCP, 275:174-196, 2014 (r-adaptação 3D rápida)
- Halila, G., Chen, G., Shi, Y. et al., Aerospace Science and Technology, vol. 91, 2019

## 8 (Parte II). Resolver com MENOS densidade, sem adaptar malha: alternativas de modelo e discretização

O ponto que faltava: adaptar malha redistribui o custo, mas não o elimina. O objetivo mais ambicioso, apontado como direção no doutoramento do SA-sLM2015 e nunca concretizado, é o modelo dar respostas corretas em malha grossa. Alternativas, da mais radical para a mais incremental:

**8.1 Mudança de variável que lineariza a frente (level-set da intermitência).** O problema não é a física, é que γ salta de 0 para 1 em poucas células. Resolve antes um campo suave φ (tipo função distância ao onset, ou o quantil da distribuição de intermitência) e recupera γ analiticamente, γ igual a g(φ) com g afiada mas exata. A malha grossa resolve φ porque φ varia devagar, e a nitidez vive na função analítica, não na malha. É a mesma jogada do level-set face ao Heaviside em interfaces. Nunca foi feito num modelo de transição. Se funcionar, o requisito de malha em corda cai para o do escoamento médio. Bónus: o gradiente adjunto de φ é suave por construção, o que ataca em simultâneo o teu problema do gradiente. Candidata a contribuição principal da tese.

**8.2 Modelar a frente em vez de a resolver (distribuição universal).** Dhawan e Narasimha (1958) mostraram que a forma de γ(x) através da frente é quase universal quando normalizada pelo onset e pelo comprimento de transição. Inverte a lógica: o transporte só determina onde começa e quanto dura, a forma é imposta analiticamente. A malha só precisa de resolver o onset, não o degrau. Próximo em espírito dos acoplamentos e^N, mas nunca embutido num LCTM local e diferenciável.

**8.3 Sub-malha embebida para as equações de transição (dual mesh).** Escoamento médio na malha grossa, γ e Reθt numa sub-malha refinada em corda alinhada à superfície, com transferência conservativa entre as duas. Análogo do wall-modeled LES e das camadas limite embebidas. As equações de transição são baratas, por isso o custo extra é pequeno, e a malha do escoamento médio volta a mandar no custo total. Complicação honesta: a diferenciação da transferência entre malhas entra na cadeia adjunta.

**8.4 Ordem alta local só no bloco de transição.** p-enriquecimento das equações de γ e Reθt na banda da frente (mantendo 2ª ordem no resto e no escoamento médio). Mesma resolução da frente com menos pontos. Meio caminho entre discretização e adaptação.

**8.5 Correção de malha grossa por ML consistente com o modelo.** Treinar um termo corretivo para o modelo em malha grossa reproduzir a frente da malha fina, usando o teu adjunto para fornecer os gradientes do treino (model-consistent training aplicado ao erro de discretização, não ao erro de modelo). A literatura FIML corrige física, ninguém corrige discretização de transição. Liga ao 6.6.

**8.6 Largura de frente como parâmetro do modelo.** Reformular os fontes para a largura da frente ser controlada pelo modelo (parâmetro físico-numérico calibrado), não um resultado acidental da malha e do esquema. Versão permanente do G3: em vez de homotopia temporária, um modelo cuja frente é sempre resolvível por construção, com calibração para recuperar os integrais certos (posição do onset e défice de fricção). É provavelmente a concretização do que o doutoramento do sLM2015 deixou em aberto.

**8.7 Convecção de 2ª ordem de γ com o controlo devolvido ao modelo.** O teu upwind de 1ª ordem espalha a frente com difusão numérica dependente da malha, incontrolável e não convergente. 2ª ordem com limitador suave, combinada com 8.6, tira a largura da frente das mãos da malha e entrega-a ao modelo. Sozinha não resolve, mas é pré-requisito de rigor para 8.1 e 8.6.

**Como isto compõe com a Parte II toda:** as alternativas de modelo (8.1, 8.2, 8.6) e a adaptação (Fase 1) não competem, multiplicam-se. O teto máximo é 8.1 mais uma r-adaptação leve: frente linearizada na variável resolvida, pontos concentrados só onde o escoamento médio precisa. E a Fase 0 (o estudo do erro do gradiente) continua a ser o alicerce de qualquer uma destas vias, porque é ela que define a métrica de sucesso: mesmo erro de gradiente com quantas vezes menos pontos.

---

# IDEIAS POR TENTAR (fila curta, mantida pela sessão de implementação)

## F1. Opt 2D com adjunto de transição congelada (frozen-transition)

Análogo do frozen-turbulence adjoint, mas congelando SÓ γ e Re̅θt na
linearização (SA continua exato): correr a otimização 2D de produção
(NACA0012) com esse gradiente aproximado e ver no que dá — onde converge,
quanto perde no ótimo vs adjunto exato, e qual o ângulo entre o gradiente
congelado e o verdadeiro (métrica do G4/5.2, que só nós conseguimos medir
por termos o adjunto exato verificado). Interesse: se o gradiente
frozen-transition for "suficientemente descendente", é um adjunto muito mais
barato e robusto (remove as linhas/colunas stiff da transição do sistema
adjunto) para as fases iniciais da opt — homotopia natural com o 5.3.
Implementação provável: máscara no dRdw/dRdwPre análoga ao `frozenTurbulence`
existente mas por-variável (7–8), + opção runtime.

## F2. Convergir em modo transiente quando o estacionário estagna

Testar se os casos SA-GR que estagnam em steady (plateaus profundos, NK
wall) convergem em modo time-accurate (unsteady/BDF do ADflow): as bolhas
de separação laminar na zona de transição podem mover-se — a solução não é
necessariamente estacionária, e nesse caso o stall do solver steady é
físico, não numérico. O que medir: (i) um caso que estagna em steady,
corrido em transiente até regime — oscila (frequência/amplitude da bolha e
da frente) ou assenta num fixed point?; (ii) se assenta: o transiente serve
de continuação para reentrar no steady/NK (warm start); (iii) se oscila: o
plateau steady corresponde à média temporal? Implicações para a opt: se o
escoamento é genuinamente instacionário, o adjunto steady está mal posto
nesses pontos — liga ao 6.5 (harmonic balance) e explica estagnações da
otimização que nenhum melhoramento de solver steady resolve.

---

# REGISTO DE IDEIAS TENTADAS (mantido pela sessão de implementação)

Cada entrada: ideia → o que se fez → conclusão concisa → onde confirmar.

## T1. G1 (parcial): reciclagem de Krylov no adjunto — LGMRES  [2026-08-22]

**O que se fez.** Opção `adjointSolver="LGMRES"` + `adjointLGMRESAugDim`
(default 2) no ADflow — LGMRES guarda `aug` vetores de erro através dos
restarts do GMRES. 4 ficheiros, 26 linhas, commit `31b467c0` deste branch;
sem tocar em equações nem AD (zero Tapenade, gradientes bit-consistentes).
Testado com o PC de produção (ADPC, ILU(3), ASM(3)) em:
- 2D NACA0012 v2 L1 (malha da opt 2D), α=2.62, adjuntos cl+cd a 1e-9,
  GMRES vs LGMRES(aug4) × subespaço {400,200,100,50} — Deucalion job 1844046;
- 3D AR5 wing+winglet L2 (1.48M células), α=5.71, primal por restart do
  volume convergido, adjuntos a 1e-8, subespaço {400,200,100} — job 1844151.

**Conclusões.**
1. Correção verificada: sensibilidades iguais entre todos os solves
   (dispersão 7e-10 rel. no 2D, 2e-9 no 3D).
2. Tempo: empate GMRES↔LGMRES em todos os subespaços que convergem (<1%).
   O ganho não é velocidade.
3. Robustez: a ss50 (2D) o GMRES entra em ciclo-limite no restart e NUNCA
   converge (2820+ it. preso a 7e-1); o LGMRES converge (1192 s). O LGMRES
   elimina o modo de falha por estagnação — seguro de custo zero.
4. Memória: real mas modesta no 3D — ss400→ss100 poupa 17–18% do RSS medido
   (3.83→3.16 GB/rank; base de Krylov ~1.2 de 3.8 GB/rank, o chão é PC+malha)
   por +38% de tempo. No 2D o RSS nem mexe (PC domina).
5. Recomendação: manter GMRES(400) onde cabe; LGMRES quando a memória
   obrigar a encolher o subespaço (c376-class: base ss400 ~350 GB, OOM job
   1826948).

**Onde confirmar.**
- Runs + PURPOSE.md + logs + result_*.json:
  `/home/mdo/Desktop/Run/MDO_PhD/Transition/gama_rethetha/02_adjoint_checks/lgmres_2d_naca0012/`
  e `.../02_adjoint_checks/lgmres_3d_ar5L2/` (espelho HPC:
  `/projects/F202500002HPCVLABISTUL/pedrocardoso/MDO_Lab_IST/TransitionModel/lgmres_{2d_naca0012,3d_ar5L2}/`).
- Slides 30–31 de `SA_gamma_Retheta_status_2026-08-22.pptx` (gerado por
  `make_deck_2026-08-22.py`).
- Código: `git show 31b467c0` neste branch (`opt-improvements-2026-08`).

**Notas colaterais úteis.**
- O padrão restart+adjunto (primal ~2 min via `restartFile`) torna A/Bs de
  adjunto 3D baratos — reutilizar.
- Armadilha de instalação: `build/lib` está tracked no repo; num clone novo
  `pip install .` empacota o pyADflow VELHO com banner novo — `rm -rf build`
  antes de instalar.


## T2. G1 (parcial): precondicionador field-split por bloco físico no adjunto  [2026-08-23]

**O que se fez.** `globalPreconditioner="field split"` +
`adjointFieldSplitType` (multiplicative default / additive) — PCFIELDSPLIT
sobre {flow ρ..ρE}, {ν̃}, {γ, Re̅θt}, cada split com ASM+ILU via options
database; matriz do PC passa a AIJ (os IS por variável cortam os blocos
8×8 da BAIJ). 5 ficheiros, commit `2a446a74`, branch `fieldsplit-adjoint-pc`
(a partir de `transition-models`); PC-only, zero Tapenade. Testado no 2D
NACA0012 v2 L1 (caso da opt 2D), jobs Deucalion 1844744+1844755, baseline
ASM do job 1844046 (T1).

**Conclusões.**
1. Correção verificada: a ss400 sensibilidades iguais ao baseline (~1e-8 rel).
2. **fsmult+LGMRES(4) @ ss400: −29% no tempo dos adjuntos cl+cd (603 vs
   851 s)**; fsmult+GMRES −19%. ~10× mais iterações KSP, mas aplicação do
   PC muito mais barata que o ILU(3)/ASM(3) monolítico inner4/outer10.
3. A ss≤200 o field-split estagna SEMPRE (mesmo com LGMRES — que no
   monolítico curava o ss50): sem o acoplamento inter-bloco no PC, o GMRES
   reiniciado não converge neste sistema. Precisa do subspace grande ⇒ não
   serve para o regime memory-constrained (c376-class); aí fica o T1.
4. additive sempre ≥ multiplicative em tempo — descartar.
5. Split a 4 {flow}{ν̃}{γ}{Re̅θt} (`adjointFieldSplitBlocks=4`, commit
   `f3116cc5`, job 1844845): empate com o split a 3 a ss400 (±2-3%) e o
   mesmo stall a ss≤200 — a escala γ↔Re̅θt não é o fator limitante; o que
   custa é cortar flow↔turbulência/transição. Default fica 3.
6. Candidato a merge no `transition-models`: ganho real onde ss400 cabe,
   opt-in, default intacto ('additive Schwarz').

**Onde confirmar.** `.../02_adjoint_checks/fieldsplit_2d_naca0012/`
(PURPOSE.md, logs, result_*.json; espelho HPC
`.../TransitionModel/fieldsplit_2d_naca0012/`). Código: `git show 2a446a74`.
