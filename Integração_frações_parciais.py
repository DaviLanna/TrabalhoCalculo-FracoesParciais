import math

# Define a precisão para a saída de decimais
PRECISAO_VISUAL = 6
TOLERANCIA = 1e-9

## === CLASSE AUXILIAR PARA POLINÔMIOS ===
class Polinomio:
    """Representa um polinômio usando uma lista de coeficientes (do maior para o menor grau)."""
    def __init__(self, coefs):
        self.coefs = self._limpar_coefs(coefs)
        self.grau = len(self.coefs) - 1

    def _limpar_coefs(self, coefs):
        inicio = 0
        while inicio < len(coefs) - 1 and coefs[inicio] == 0:
            inicio += 1
        return [self._arredondar_para_zero(c) for c in coefs[inicio:]]

    def _arredondar_para_zero(self, valor, tolerancia=TOLERANCIA):
        # Arredonda valores muito próximos de zero
        if abs(valor) < tolerancia: return 0.0
        return valor
    
    def _formata_coeficiente_decimal(self, valor, precisao=PRECISAO_VISUAL):
        # Formata para decimal com precisão definida
        if abs(valor) == 1.0: return ""
        return f"{valor:,.{precisao}f}".rstrip('0').rstrip('.')

    def __str__(self):
        """Representação de string formatada para o polinômio."""
        termos = []
        n = self.grau
        for i, coef in enumerate(self.coefs):
            if coef == 0:
                continue
            
            p = n - i
            
            # Coeficiente
            abs_coef = abs(coef)
            sinal = "+" if coef > 0 and i > 0 else "-" if coef < 0 else ""
            
            coef_str = ""
            if p == 0: # Termo constante
                 coef_str = f"{abs_coef:g}"
            elif abs_coef == 1.0: # Coeficiente unitário
                 coef_str = ""
            else:
                 coef_str = f"{abs_coef:g}"
                 
            # Variável x
            if p == 1:
                var_str = "x"
            elif p > 1:
                var_str = f"x^{p}"
            elif p == 0:
                var_str = ""
            
            termos.append(f"{sinal}{coef_str}{var_str}")
        
        # Ajusta o sinal do primeiro termo
        if not termos: return "0"
        primeiro_termo = termos[0]
        if primeiro_termo.startswith('+'): termos[0] = primeiro_termo[1:]

        return "".join(termos) or "0"

## === FUNÇÕES DE DECOMPOSIÇÃO EM FRAÇÕES PARCIAIS (FP) ===

def resolver_sistema_linear(A, b):
    """
    Resolve um sistema de equações lineares Ax = b usando Eliminação Gaussiana.
    """
    n = len(A)
    M = [A[i] + [b[i]] for i in range(n)] 

    for i in range(n):
        linha_max = i
        for k in range(i + 1, n):
            if abs(M[k][i]) > abs(M[linha_max][i]): linha_max = k
        M[i], M[linha_max] = M[linha_max], M[i]

        if M[i][i] == 0: raise ValueError("Sistema singular detectado.")

        for k in range(i + 1, n):
            fator = M[k][i] / M[i][i]
            for j in range(i, n + 1): M[k][j] -= fator * M[i][j]

    x = [0] * n
    for i in range(n - 1, -1, -1):
        soma = sum(M[i][j] * x[j] for j in range(i + 1, n))
        x[i] = (M[i][n] - soma) / M[i][i]

    pol = Polinomio(x)
    return pol.coefs

def decomposicao_fracoes_parciais(P, Q, fatores):
    """Realiza a decomposição focada nos tipos de problemas propostos."""
    if P.grau >= Q.grau: raise NotImplementedError("Não suporta frações impróprias.")
    
    if Q.grau == 2 and fatores[0][0] == 'linear' and fatores[1][0] == 'linear':
        x1, x2 = fatores[0][1][0], fatores[1][1][0]
        A_matriz = [ [1.0, 1.0], [-x2, -x1] ]
        b_vetor = P.coefs
        constantes = resolver_sistema_linear(A_matriz, b_vetor)
        A_const, B_const = constantes
        return [
            ('log', A_const, Polinomio([1, -x1])),
            ('log', B_const, Polinomio([1, -x2])),
        ]
        
    elif Q.grau == 4 and fatores[0][0] == 'quadratic_irr' and fatores[1][0] == 'quadratic_irr':
        c1 = fatores[0][1][2]
        c2 = fatores[1][1][2]
        A_matriz = [
            [1.0, 0.0, 1.0, 0.0], [0.0, 1.0, 0.0, 1.0], 
            [c2, 0.0, c1, 0.0], [0.0, c2, 0.0, c1] 
        ]
        b_vetor = [0.0, 0.0] + P.coefs
        constantes = resolver_sistema_linear(A_matriz, b_vetor)
        A_const, B_const, C_const, D_const = constantes
        return [
            ('quadratic', (A_const, B_const), Polinomio([1, 0, c1])),
            ('quadratic', (C_const, D_const), Polinomio([1, 0, c2])),
        ]

    elif Q.grau == 2 and len(fatores) == 1 and fatores[0][0] == 'quadratic_irr':
        return [
            ('quadratic_direta', tuple(P.coefs), Q)
        ]
    
    raise NotImplementedError("Estrutura do denominador não suportada.")
    
def integrar_fracao(tipo, numerador, denominador):
    """
    Integra uma única fração parcial e aplica formatação limpa (decimais).
    """
    pol_formatter = Polinomio([])
    
    # Caso 1: Logarítmica $\int A/(x-x_1)dx$
    if tipo == 'log':
        A = numerador
        if A == 0: return "0"
        
        A_str = pol_formatter._formata_coeficiente_decimal(abs(A))
        sinal = "-" if A < 0 else ""
        
        return f"{sinal}{A_str}ln(|{denominador}|)"

    # Caso 2/4: Quadrática/Arctan $\int (Ax+B)/(ax^2+bx+c)dx$
    elif tipo == 'quadratic' or tipo == 'quadratic_direta':
        A, B = numerador
        a, b_q, c = denominador.coefs
        
        # Termo 1 (ln)
        coef_ln = pol_formatter._arredondar_para_zero(A / (2 * a))
        C_const = pol_formatter._arredondar_para_zero(B - coef_ln * b_q)
        
        resultado = ""
        
        if coef_ln != 0:
            coef_ln_str = pol_formatter._formata_coeficiente_decimal(abs(coef_ln))
            sinal = "-" if coef_ln < 0 else ""
            resultado += f"{sinal}{coef_ln_str}ln({denominador})"
            
        # Termo 2 (arctan)
        if C_const != 0:
            h = pol_formatter._arredondar_para_zero(b_q / (2 * a))
            k_quadrado = c / a - h**2
            
            if k_quadrado <= 0: raise ValueError("O fator quadrático não é irredutível.")
                
            k = math.sqrt(k_quadrado)
            coef_arctan = pol_formatter._arredondar_para_zero(C_const / (a * k))
            
            if coef_arctan != 0:
                sinal_sep = " - " if coef_arctan < 0 and resultado else " + " if coef_arctan > 0 and resultado else "-" if coef_arctan < 0 else ""
                
                abs_coef_arctan = abs(coef_arctan)
                coef_arctan_str = pol_formatter._formata_coeficiente_decimal(abs_coef_arctan)
                
                # Saída limpa: usa decimais nos argumentos e evita x+0 ou /1
                k_str = f"{k:,.{PRECISAO_VISUAL}f}".rstrip('0').rstrip('.')
                x_parte = f"x+{h:g}" if h != 0 else "x"
                
                termo_arctan = f"{sinal_sep}{coef_arctan_str}arctan({x_parte}/{k_str})"
                resultado += termo_arctan
        
        return resultado.strip(" +").strip(" -")


## === FUNÇÃO PRINCIPAL DE INTEGRAÇÃO ===

def integrar_funcao_racional(coefs_P, coefs_Q, teste_num):
    """
    Função principal para integrar P(x)/Q(x) usando FP, com saída formatada para terminal simples.
    """
    P = Polinomio(coefs_P)
    Q = Polinomio(coefs_Q)
    
    # --- SAÍDA DO CABEÇALHO ---
    print(f"\n============================================================")
    print(f"## 🧪 Teste {teste_num}: Integração de P(x)/Q(x)")
    print(f"Função Racional: ({P}) / ({Q})")
    print(f"--- Análise e Decomposição ({P.grau} / {Q.grau}) ---")

    # 1. Identificar Fatores
    fatores = []
    
    if Q.grau == 2 and Q.coefs[0] == 1 and Q.coefs[2] < 0:
        b_q, c_q = Q.coefs[1], Q.coefs[2]
        delta = b_q**2 - 4*c_q
        if delta > 0:
            x1 = (-b_q + math.sqrt(delta)) / 2
            x2 = (-b_q - math.sqrt(delta)) / 2
            fatores.append(('linear', (x1,)))
            fatores.append(('linear', (x2,)))
            print(f"* Fatores: (x - {P._arredondar_para_zero(x1):g})(x + {P._arredondar_para_zero(abs(x2)):g})")

    elif Q.grau == 2 and (Q.coefs[1]**2 - 4*Q.coefs[0]*Q.coefs[2]) < 0:
        fatores.append(('quadratic_irr', tuple(Q.coefs)))
        print(f"* Fator: {Q} (Irredutível)")

    elif Q.grau == 4 and Q.coefs == [1, 0, Q.coefs[2], 0, Q.coefs[4]] and Q.coefs[2] > 0 and Q.coefs[4] > 0:
        c1_mais_c2 = Q.coefs[2]
        c1_vezes_c2 = Q.coefs[4]
        delta = c1_mais_c2**2 - 4*c1_vezes_c2
        if delta >= 0:
            c1 = (c1_mais_c2 + math.sqrt(delta)) / 2
            c2 = (c1_mais_c2 - math.sqrt(delta)) / 2
            fatores.append(('quadratic_irr', (1, 0, c1)))
            fatores.append(('quadratic_irr', (1, 0, c2)))
            print(f"* Fatores: (x^2 + {P._arredondar_para_zero(c1):g})(x^2 + {P._arredondar_para_zero(c2):g})")

    if not fatores: raise ValueError("Não foi possível identificar/fatorar Q(x).")

    # 2. Decomposição
    fracoes = decomposicao_fracoes_parciais(P, Q, fatores)
    
    decomp_str = []
    for tipo, num, den in fracoes:
        if tipo == 'log':
            num_str = f"{num:,.{PRECISAO_VISUAL}f}".rstrip('0').rstrip('.')
            decomp_str.append(f"({num_str}) / ({den})")
        elif tipo == 'quadratic' or tipo == 'quadratic_direta':
            num_A_str = f"{num[0]:g}" if abs(num[0]) != 1 else ""
            num_B_str = f"{num[1]:,.{PRECISAO_VISUAL}f}".rstrip('0').rstrip('.') if num[1] != 0 else ""
            
            termo_x = f"{num[0]:g}x" if num[0] != 0 else ""
            termo_const = f"{'+' if num[1]>0 and num[0]!=0 else ''}{num_B_str}" if num[1] != 0 else ""
            num_final = f"{termo_x}{termo_const}".replace("+-", " - ")
            if num_final.startswith('+'): num_final = num_final[1:]

            decomp_str.append(f"({num_final}) / ({den})")
            
    print(f"--- 📝 Decomposição em Frações Parciais ---")
    print(f"({P}) / ({Q}) = {' + '.join(decomp_str)}")


    # 3. Integração Termo a Termo
    print(f"--- ✨ Integração Termo a Termo ---")
    integral_results = []
    
    for i, (tipo, num, den) in enumerate(fracoes):
        integral_termo = integrar_fracao(tipo, num, den)
        integral_results.append(integral_termo)
        
        if tipo == 'log':
             frac_str = f"({num:g}) / ({den})"
        elif tipo == 'quadratic' or tipo == 'quadratic_direta':
             frac_str = f"({num[0]:g}x+{num[1]:g}) / ({den})"
             
        # Saída simplificada para o termo
        print(f"* Integral de {frac_str} = {integral_termo}")
            
    # 4. Resultado Final
    print(f"--- ✅ Resultado Final ---")
    final_result = "".join(f" {('+' if res[0] != '-' and i > 0 else '')} {res}" for i, res in enumerate(integral_results)).strip()
    
    print(f"Resultado Final: {final_result} + C")
    
    print("\n" + "="*50)
    return final_result


## === TESTES (Exemplos Propostos no Trabalho) ===

def rodar_exemplos():
    """Executa os exemplos para validação."""
    
    # Exemplo 1: $\int\frac{3x+5}{x^{2}-4}dx$
    integrar_funcao_racional([3, 5], [1, 0, -4], 1)

    # Exemplo 2: $\int\frac{x-1}{x^{2}+2x+10}dx$
    integrar_funcao_racional([1, -1], [1, 2, 10], 2)

    # Exemplo 3: $\int\frac{x+1}{x^4+3x^2+2}dx$
    integrar_funcao_racional([1, 1], [1, 0, 3, 0, 2], 3)


if __name__ == "__main__":
    rodar_exemplos()