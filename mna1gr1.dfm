object Form1: TForm1
  Left = 177
  Top = 73
  Width = 1016
  Height = 559
  Caption = 'RespFreq_CMOS'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  Menu = MainMenu1
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Memo1: TMemo
    Left = 0
    Top = 0
    Width = 1000
    Height = 500
    Align = alClient
    Font.Charset = ANSI_CHARSET
    Font.Color = clWindowText
    Font.Height = -11
    Font.Name = 'Courier New'
    Font.Style = []
    Lines.Strings = (
      'Programa do trabalho de 2016.1'
      ''
      
        'An'#225'lise de ponto de opera'#231#227'o e de resposta em frequ'#234'ncia de circ' +
        'uitos lineares contendo transistores MOS'
      
        'Baseado no programa "mna1" do professor Ant'#244'nio Carlos M. de Que' +
        'iroz - acmq@ufrj.br'
      ''
      'Por: Bruno Granato'
      '     Jo'#227'o Felipe Guedes'
      '     Lucas de Andrade Cerqueira')
    ParentFont = False
    ScrollBars = ssBoth
    TabOrder = 0
  end
  object OpenDialog1: TOpenDialog
    FileName = 
      'C:\Users\Lucas Cerqueira\Documents\EngEletronica\5o periodo\Circ' +
      'uitos Eletricos 2\Trabalho\Circuitos\ampcmos.net'
    Filter = 'netlists|*.net'
    Options = [ofHideReadOnly, ofNoChangeDir, ofEnableSizing]
    Left = 632
    Top = 16
  end
  object MainMenu1: TMainMenu
    Left = 88
    Top = 16
    object Arquivo1: TMenuItem
      Caption = 'Arquivo'
      object Abrir1: TMenuItem
        Caption = 'Abrir'
        OnClick = Abrir1Click
      end
    end
    object Opes1: TMenuItem
      Caption = 'Op'#231#245'es'
      object MostrarEstampasDC1: TMenuItem
        AutoCheck = True
        Caption = 'Exibir estampas DC'
      end
      object MostrarEstampasAC1: TMenuItem
        Caption = 'Exibir estampas AC'
      end
      object MostrarParametrosMOS1: TMenuItem
        AutoCheck = True
        Caption = 'Exibir par'#226'metros MOS'
      end
      object IniciarMOSconduzindo1: TMenuItem
        AutoCheck = True
        Caption = 'Iniciar MOS conduzindo'
      end
    end
  end
end
