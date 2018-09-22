
import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
import datetime
import time

def main():

    sair = False
    while(not sair):
        
        menu = int(input("\n\nSelecione vizualizar testes de transformada(1) Testes de dados(2) ou Resolução do problema proposto (3): "))
        
        if menu == 1:
            
            menu_teste = int(input("\n\nSelecione teste a) (1), b) (2) ou c) (3): "))


            if menu_teste == 1:

                print("\n\nTeste a) f = [5,-1,3,1]",end="\n\n")
                f = [5,-1,3,1]

                print("Implementação ineficiente da Transformada de Fourier direta",end="\n")
                c1 = fft_direta(f)
                print_vector(c1)

                print("Implementação recursiva da FFT direta",end="\n")
                c = [0,0,0,0]
                fftrec(c,f,2,True)
                print_vector(np.divide(c,4))

                print("\n\nFazendo c = (2, 1/2 + i/2, 2, 1/2 − i/2)",end="\n\n")
                c = [complex(2,0),complex(0.5,0.5),complex(2,0),complex(0.5,-0.5)]
                
                print("Implementação ineficiente da Transformada de Fourier inversa",end="\n")
                f1 = fft_inversa(c)
                print_vector(f1)

                print("Implementação recursiva da FFT inversa",end="\n")
                f = [0,0,0,0]
                fftrec(f,c,2,False)
                print_vector(f)
                
            if menu_teste == 2:

                print("\n\nTeste b) f = [6, 2, 5, 2, 11, 2, 8, 8]",end="\n\n")
                f = [6, 2, 5, 2, 11, 2, 8, 8]

                print("Implementação ineficiente da Transformada de Fourier direta",end="\n")
                c1 = fft_direta(f)
                print_vector(c1)

                print("Implementação recursiva da FFT direta",end="\n")
                c = [0,0,0,0,0,0,0,0]
                fftrec(c,f,4,True)
                print_vector(np.divide(c,8))

                print("Implementação recursiva da FFT direta pelo numpy (rfft)",end="\n")
                print("Note que são mostrados apenas os coeficientes de frequencias positivas pois o input é real")
                c2 = np.fft.rfft(f)
                print_vector(np.divide(c2,8))
                
                print("Implementação ineficiente da Transformada de Fourier inversa",end="\n")
                f1 = fft_inversa(c1)
                print_vector(f1)

                print("Implementação recursiva da FFT inversa",end="\n")
                f3 = [0,0,0,0,0,0,0,0]
                fftrec(f3,np.divide(c,8),4,False)
                print_vector(f3)

                print("Implementação recursiva da FFT inversa pelo numpy (rfft)",end="\n")
                f2 = np.fft.irfft(c2)
                print_vector(f2)

            if menu_teste == 3:

                print("\n\nTeste b) F(x) = 10sin(x) + 7cos(30x) + 11sin(352x) − 8cos(711x) com xj = jπ/512, j = 0, 1023",end="\n\n")
                dois_n = 1024
                f = [ 0 for j in range(dois_n)]
                x = [j*math.pi/(dois_n/2) for j in range(dois_n)]
                f = [10*math.sin(x[j])+7*math.cos(30*x[j])+11*math.sin(352*x[j])-8*math.cos(711*x[j]) for j in range(dois_n)]

                
                print("Implementação ineficiente da Transformada de Fourier direta",end="\n")
                c1 = fft_direta(f)

                print("Implementação recursiva da FFT direta",end="\n")
                c = [0 for i in range(dois_n)]
                fftrec(c,f,int(dois_n/2),True)

                print("Implementação recursiva da FFT direta pelo numpy (rfft)",end="\n")
                print("Note que são mostrados apenas os coeficientes de frequencias positivas pois o input é real")
                c2 = np.fft.rfft(f)
                
                print("Implementação ineficiente da Transformada de Fourier inversa",end="\n")
                f1 = fft_inversa(c1)

                print("Implementação recursiva da FFT inversa",end="\n")
                f3 = [0 for i in range(dois_n)]
                fftrec(f3,np.divide(c,dois_n),int(dois_n/2),False)

                print("Implementação recursiva da FFT inversa pelo numpy (rfft)",end="\n")
                f2 = np.fft.irfft(c2)


                plt.figure(1)
                ax1 = plt.subplot(231)
                t = [j for j in range(int(dois_n))]
                plt.plot(t, np.absolute(c1), 'bo')
                ax1.set_title('TF Ineficiente')
                k = [j for j in range(dois_n)]
                plt.subplot(234)
                plt.plot(k, f1, 'r-')

                ax2 = plt.subplot(232)
                t = [j for j in range(int(dois_n))]
                plt.plot(t, np.absolute(np.divide(c,dois_n)), 'bo')
                ax2.set_title('FFT recursiva')
                k = [j for j in range(dois_n)]
                plt.subplot(235)
                plt.plot(k, f3, 'r-')

                ax3 = plt.subplot(233)
                t = [j for j in range(int(dois_n/2+1))]
                plt.plot(t, np.absolute(np.divide(c2,dois_n)), 'bo')
                ax3.set_title('FFT numpy')
                k = [j for j in range(dois_n)]
                plt.subplot(236)
                plt.plot(k, f2, 'r-')

                plt.show()


        if menu == 2:
            
            Alum = read_file("Alum.csv")
            nome = 'Aluminio'
            print("===============================")
            tipo_tf = int(input("Selecione transformada Lenta N2(0), recursiva implementada(1) ou numpy(2): "))
            tipo_filtro = int(input("Selecione filtro passa baixas(0) ou passa altas(1): "))
            K = int(input("Digite o indice K de corte do filtro: "))
            print("===============================")
            
            Alum_preco = [x[2] for x in Alum]

            if tipo_tf == 0:
                timedir,timeinv = plotComTF_Lenta(Alum_preco, nome, K, tipo_filtro)
            if tipo_tf == 1:
                timedir,timeinv = plotComTF_rec(Alum_preco, nome, K, tipo_filtro)
            if tipo_tf == 2:
                timedir,timeinv = plotComTF_numpy(Alum_preco, nome, K, tipo_filtro)

            print("===============================")
            print("\n\nTempo decorrido para transformada direta: %f ms" %(timedir))
            print("\nTempo decorrido para transformada inversa: %f ms" %(timeinv))   
            print("\nDesvio Padrao da serie temporal: %f" %(desvPad(Alum)))
            print("===============================")


        if menu == 3:

            Cattle = read_file("Cattle.csv")

            for i in range(30):
                print(Cattle[i])
                
                        
        sair = bool(int(input("Sair: 1 Continuar: 0\n")))

    

def fft_inversa(c):
    n = int(len(c)/2)
    f = [0 for i in range(2*n)]
    
    for j in range(0,2*n):

        for k in range(0,2*n):

            f[j] += c[k]*cmath.exp(complex(0,k*j*math.pi/n))

    return f

def fft_direta(f):
    n = int(len(f)/2)
    c = [0 for i in range(2*n)]
    
    for k in range(0,2*n):

        for j in range(0,2*n):

            c[k] += (f[j]*cmath.exp(complex(0,-k*j*math.pi/n)))/(2*n)

    return c
        
def fftrec(c,f,n,direto):
    fe = [0 for i in range(2*n)]
    fo = [0 for i in range(2*n)]
    even = [c[2*j] for j in range(n)]
    odd = [c[2*j+1] for j in range(n)]
    
    if n == 1:

        c[0] = f[0]+f[1]
        c[1] = f[0]-f[1]

    else:

        for j in range(n):

            fe[j] = f[2*j]
            fo[j] = f[2*j+1]

        fftrec(even,fe,int(n/2),direto)
        fftrec(odd,fo,int(n/2),direto)

        for j in range(n):

            if direto:
                eij = cmath.exp(complex(0,-j*math.pi/n))

            else:
                eij = cmath.exp(complex(0,j*math.pi/n))

            c[j] = even[j] + eij*odd[j]
            c[j+n] = even[j] - eij*odd[j]

#=======================================================================

def print_vector(vetor):
    ''' matriz -> null '''
    for k in range(45):
        print("=", end="")
    print("",end="\n")
    for n in range (0,len(vetor),1):
        print ("| (%.5f, %.5fi)" %(vetor[n].real, vetor[n].imag), end="")
    print("|")
    for k in range(45):
        print("=", end="")
    print("", end="\n\n")
    return

#=======================================================================

def read_file(file):
    '''
        Lê arquivo
    '''
    f = open(file)
    arquivo = [x.strip('\n') for x in f.readlines()]
    arquivo = [x.split(', ') for x in arquivo]

    for l in range(len(arquivo)):
        for i in range(3):
            if i == 0:
                arquivo[l][i] = datetime.date(int(arquivo[l][i][0:4]),
                                              int(arquivo[l][i][4:6]),
                                              int(arquivo[l][i][6:9]))
            if i == 1:
                arquivo[l][i] = int(arquivo[l][i])

            if i == 2:
                arquivo[l][i] = float(arquivo[l][i])
                
    return arquivo

#=======================================================================

def mean(numbers):
    '''
        Calcula media de vetor
    '''
    return float(sum(numbers)) / max(len(numbers), 1)

#=======================================================================

def desvPad(file):
    '''
        Calcula desvio padrao
    '''

    valores = [x[2] for x in file]
    media = mean(valores)
    argsum = np.divide(np.power(np.add(valores,-media),2),len(valores))
    desvPad = math.sqrt( sum( argsum ) )
                
    return float(desvPad)

#=======================================================================

def passa_baixa(c,K):
    c2 = [c[j] if j < K-1 else 0 for j in range(len(c))]
    return c2

def passa_alta(c,K):
    c2 = [c[j] if j >= K-1 else 0 for j in range(len(c))]
    return c2

def passa_baixa2(c,K):
    tam = len(c)
    c2 = [0 for j in range(tam)]
    for j in range(tam):
        if (j < K-1) or (j > tam-K-1 ):
            c2[j] = c[j]
    return c2

def passa_alta2(c,K):
    tam = len(c)
    c2 = [0 for j in range(tam)]
    for j in range(tam):
        if (j >= K-1) and (j <= tam-K-1 ):
            c2[j] = c[j]
    return c2

#=======================================================================

def plotComTF_Lenta(valor, nome, K, tipo_filtro):
    '''
        valor: vetor de floats
        nome: string
        K: inteiro
        tipo_filtro: inteiro
    '''
    print("\n\nNote que a frequencia de indice 0 está omitida do plot")
    time_dir_st = time.time()
    c = fft_direta(valor)
    time_dir_end = time.time()
    time_dir = (time_dir_end - time_dir_st)*1000
    
    plt.figure(1)
    ax1 = plt.subplot(211)
    t = [j for j in range(1,len(valor))]
    plt.plot(t, np.absolute(c)[1:len(c)], 'bo')
    ax1.set_title(nome)
    plt.subplot(212)
    k = [j for j in range(len(valor))]
    plt.plot(k, valor, 'b-')

    if tipo_filtro == 0:
        c_filtrado = passa_baixa2(c,K)

    if tipo_filtro == 1:
        c_filtrado = passa_alta2(c,K)
    
    time_inv_st = time.time()
    f_filtrado = fft_inversa(c_filtrado)
    time_inv_end = time.time()
    time_inv = (time_inv_end - time_inv_st)*1000

    
    
    plt.subplot(212)
    plt.plot(k, np.real(f_filtrado), 'r-')
    plt.show()

    return (time_dir,time_inv)

#=======================================================================

def plotComTF_rec(valor, nome, K, tipo_filtro):
    '''
        valor: vetor de floats
        nome: string
        K: inteiro
        tipo_filtro: inteiro
    '''
    print("\n\nNote que a frequencia de indice 0 está omitida do plot")
    c = [0 for i in range(len(valor))]
    time_dir_st = time.time()
    fftrec(c,valor,int(len(valor)/2),True)
    time_dir_end = time.time()
    time_dir = (time_dir_end - time_dir_st)*1000
    
    plt.figure(1)
    ax1 = plt.subplot(211)
    t = [j for j in range(1,len(valor))]
    plt.plot(t, np.absolute(c)[1:len(c)], 'bo')
    ax1.set_title(nome)
    plt.subplot(212)
    k = [j for j in range(len(valor))]
    plt.plot(k, valor, 'b-')

    if tipo_filtro == 0:
        c_filtrado = passa_baixa2(c,K)

    if tipo_filtro == 1:
        c_filtrado = passa_alta2(c,K)
    
    valor_filtrado = [0 for i in range(len(valor))]
    time_inv_st = time.time()
    fftrec(valor_filtrado,np.divide(c_filtrado,len(c_filtrado)),int(len(valor)/2),False)
    time_inv_end = time.time()
    time_inv = (time_inv_end - time_inv_st)*1000

    
    
    plt.subplot(212)
    plt.plot(k, np.real(valor_filtrado), 'r-')
    plt.show()

    return (time_dir,time_inv)

#=======================================================================

def plotComTF_numpy(valor, nome, K, tipo_filtro):
    '''
        valor: vetor de floats
        nome: string
        K: inteiro
        tipo_filtro: inteiro
    '''
    print("\n\nNote que a frequencia de indice 0 está omitida do plot")
    time_dir_st = time.time()
    c = np.fft.rfft(valor)
    time_dir_end = time.time()
    time_dir = (time_dir_end - time_dir_st)*1000
    
    plt.figure(1)
    ax1 = plt.subplot(211)
    t = [j for j in range(1,int(len(valor)/2)+1)]
    plt.plot(t, np.absolute(c)[1:], 'bo')
    ax1.set_title(nome)
    k = [j for j in range(len(valor))]
    plt.subplot(212)
    plt.plot(k, valor, 'b-')
    
    if tipo_filtro == 0:
        c_filtrado = passa_baixa(c,K)

    if tipo_filtro == 1:
        c_filtrado = passa_alta(c,K)

    time_inv_st = time.time()
    valor_filtrado = np.fft.irfft(c_filtrado)
    time_inv_end = time.time()
    time_inv = (time_inv_end - time_inv_st)*1000
    
    plt.subplot(212)
    plt.plot(k, valor_filtrado, 'r-')
    plt.show()

    return (time_dir,time_inv)

    
