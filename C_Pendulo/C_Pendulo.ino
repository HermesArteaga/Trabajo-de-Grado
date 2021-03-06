//////////////////////////
//Librerias del programa//
//////////////////////////

#include <SPI.h>
#include <Adafruit_GFX.h>
#include <Adafruit_PCD8544.h>
Adafruit_PCD8544 display = Adafruit_PCD8544(52, 50, 48, 46, 44);
//Adafruit_PCD8544 display = Adafruit_PCD8544(44, 46, 48, 50, 52);

//////////////////////////
//Variables del programa//
//////////////////////////

volatile unsigned int conteo_b = 32768;
volatile unsigned int conteo_p = 32768;
int pulsos_b=0;
int pulsos_p=0;
int pwm_valor=0;
byte  estado=0;
int  muestras=3;
byte puntero_p=0;
byte puntero_b=0;
float grados_b=0.0;
float grados_p=0.0;
float arreglo_p[5]={0.0, 0.0, 0.0, 0.0, 0.0};
float arreglo_b[5]={0.0, 0.0, 0.0, 0.0, 0.0};
float media_p=0.0;
float media_b=0.0;
float dato;

///////////////////////
//Definicion de Pines//
///////////////////////

#define inicio  22
#define reset   24
#define derecha  26
#define izquierda  28
#define led_derecha  30
#define led_izquierda  32
#define motor_a  9
#define motor_b  8
#define ena_motor  10
#define pwm_motor  11

//////////////////////////
//Configuracion de pines//
//////////////////////////
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
void setup() {
  // put your setup code here, to run once:
  //TCCR1B = TCCR1B & 0b000 | 0x02;
  
  Serial.begin(115200);  

  pinMode(inicio,INPUT);
  pinMode(reset,INPUT);
  pinMode(derecha,INPUT);
  pinMode(izquierda,INPUT);
  pinMode(led_derecha,OUTPUT);
  pinMode(led_izquierda,OUTPUT);
  pinMode(motor_a,OUTPUT);
  pinMode(motor_b,OUTPUT);
  pinMode(ena_motor,OUTPUT);
  pinMode(pwm_motor,OUTPUT);

  digitalWrite(motor_a,LOW);
  digitalWrite(motor_b,LOW);
  digitalWrite(ena_motor,LOW);
  analogWrite(pwm_motor,pwm_valor);

  display.begin();
  display.setContrast(60);
  display.clearDisplay();
  display.display();
  //display.fillRect(0, 0, display.width(), display.height(), 1);
  //display.setTextColor(WHITE, BLACK); // 'inverted' text
  display.setTextSize(1);
  display.setCursor(9,0);
  display.print("Universidad");
  display.setCursor(15,8);
  display.print("de Nari");
  display.write(164);
  display.write('o');
  //display.setCursor(12,24);
  //display.setTextSize(2);
  //display.println("GIIEE");  
  display.setCursor(0,24);
  display.println("Villota Brian");
  display.println("Arteaga Hermes");
  display.display();
  delay(2000);
  display.clearDisplay();
  display.print("  Pendulo de      Furuta");
  display.setCursor(24,24);
  display.print("Grados");
  display.setCursor(0,32);
  display.print("Pendulo=");
  display.setCursor(0,40);
  display.print("Brazo=");
  display.display();

  pwm_valor=40;
    
  pinMode(3,INPUT);
  pinMode(2,INPUT);
  pinMode(18,INPUT);
  pinMode(19,INPUT);
  attachInterrupt(digitalPinToInterrupt(18), eb_A, RISING);  //Bajo a Alto
  attachInterrupt(digitalPinToInterrupt(19), eb_B, RISING);  //
  attachInterrupt(digitalPinToInterrupt(2), ep_A, RISING);  //
  attachInterrupt(digitalPinToInterrupt(3), ep_B, RISING);  //
  
}

//////////////////////
//Programa Principal//
//////////////////////

void loop() {
  // put your main code here, to run repeatedly:
  switch (estado)
  {
    //Resposo
    case 0:
      pulsos_b=32768-conteo_b;
      grados_b=0.375*pulsos_b;
      pulsos_p=32768-conteo_p;
      grados_p=0.12*pulsos_p;
      men_grados();
      
      if(digitalRead(derecha)==1)
      {
        estado=2;
        delay(20);
        digitalWrite(motor_b,HIGH);
        digitalWrite(motor_a,LOW);
        digitalWrite(ena_motor,HIGH);
        analogWrite(pwm_motor,pwm_valor);
        digitalWrite(led_derecha,HIGH);
      }else if(digitalRead(izquierda)==1)
      {
        estado=3;
        delay(20);
        digitalWrite(motor_a,HIGH);
        digitalWrite(motor_b,LOW);
        digitalWrite(ena_motor,HIGH);
        analogWrite(pwm_motor,pwm_valor);
        digitalWrite(led_izquierda,HIGH);        
      }else if(digitalRead(reset)==1)
      {
        estado=4;
        detener();
        
      }else if (digitalRead(inicio)==1)
      {
        delay(1000);
        analogWrite(pwm_motor,0);
        estado=1;
        display.fillRect(0, 24, 84, 40, 0);
        display.setCursor(0,32);
        display.print("Enviando datos");
        display.setCursor(12,40);
        display.print("a Simulink");
        display.display();
        
        digitalWrite(motor_b,LOW);
              //delay(2);
        digitalWrite(motor_a,HIGH);
              
        digitalWrite(ena_motor,HIGH);
            //pwm_valor = int(dato);
        analogWrite(pwm_motor,50);
        delay(15);
          
        analogWrite(pwm_motor,0);
        
        
      }else if(Serial.available() > 0)
      {
        
        String str = Serial.readStringUntil('\n');
        dato = str.toFloat();
        
        //display.setCursor(0,24);
        //display.print(dato,DEC);
        //display.display();
        
        
        //if (dato=='a')
        {
          for (int i=0;i<muestras;i++)
          {
            arreglo_p[i]=0.0;
            arreglo_b[i]=0.0;
          }
          puntero_b=0;
          puntero_p=0;
          delay(2500);
          Serial.flush();
          conteo_b = 32768;
          pulsos_b=32768-conteo_b;
          grados_b=0.375*pulsos_b;
          conteo_p = 32768;
          pulsos_p=32768-conteo_p;
          grados_p=0.12*pulsos_p;

          analogWrite(pwm_motor,0);
          estado=1;
          display.fillRect(0, 24, 84, 40, 0);
          display.setCursor(0,32);
          display.print("Enviando datos");
          display.setCursor(12,40);
          display.print("a Simulink");
          display.display();
          
          digitalWrite(motor_b,LOW);
              //delay(2);
          digitalWrite(motor_a,HIGH);
              
          digitalWrite(ena_motor,HIGH);
              //pwm_valor = int(dato);
          analogWrite(pwm_motor,50);
          delay(15);
          
          analogWrite(pwm_motor,0);
          
        }
      }
      
      //delay(100);
    break;
    
    //Envio datos
    case 1:
      
      pulsos_b=32768-conteo_b;
      grados_b=0.375*pulsos_b;
      
      pulsos_p=32768-conteo_p;
      grados_p=0.12*pulsos_p;
      
      //Media brazo
      if(puntero_b<(muestras-1))
      {
        arreglo_b[puntero_b]=grados_b;
        puntero_b++;
        //Serial.print(arreglo_b[puntero_b]);
        //Serial.print(',');
      }else
      {
        arreglo_b[puntero_b]=grados_b;
        puntero_b=0;
        //Serial.println(arreglo_b[puntero_b]);
      }
      media_b=0;
      for(int i=0;i<muestras;i++)
      {
        media_b=media_b+arreglo_b[i];
      }
      media_b=media_b/muestras;
      //Serial.println(media_b);
      
      //Media pendulo
      if(puntero_p<(muestras-1))
      {
        arreglo_p[puntero_p]=grados_p;
        puntero_p++;
      }else
      {
        arreglo_p[puntero_p]=grados_p;
        puntero_p=0;
      }
      media_p=0;
      for(int i=0;i<muestras;i++)
      {
        media_p=media_p+arreglo_p[i];
      }
      media_p=media_p/muestras;
      
      //Envio informacion
      Serial.print(media_b);
      Serial.print(',');
      Serial.println(grados_p);
      
      if(Serial.available() > 0)
      {
        
        String str = Serial.readStringUntil('\n');
        dato = str.toFloat();
        if(dato<0)
        {
            if(dato<-12)
            {
              dato=12;
            }
            //digitalWrite(motor_a,HIGH);
            dato=dato*(-1);
            pwm_valor=map(dato,0,12,14,255);
            //pwm_valor=map(dato,0,12,0,255);
            digitalWrite(motor_b,HIGH);
            //delay(2);
            digitalWrite(motor_a,LOW);
            
            digitalWrite(ena_motor,HIGH);
            //pwm_valor = int(dato);
            analogWrite(pwm_motor,pwm_valor);
        }else if (dato>0)
        {
            //digitalWrite(motor_b,HIGH);
            
            if(dato>12)
            {
              dato=12;
            }
            dato=dato*1;
            pwm_valor=map(dato,0,12,14,255);
            //|||||pwm_valor=map(dato,0,12,0,255);
            digitalWrite(motor_a,HIGH);
            //delay(2);
            digitalWrite(motor_b,LOW);
            
            digitalWrite(ena_motor,HIGH);
            //pwm_valor = int(dato);
            analogWrite(pwm_motor,pwm_valor);
        }else if(dato==0)
        {
          digitalWrite(motor_a,LOW);
            //delay(2);
          digitalWrite(motor_b,LOW);
          analogWrite(pwm_motor,0);
        }
        
      }
      if(digitalRead(reset)==1)
      {
        estado=0;
        display.fillRect(0, 24, 84, 40, 0);
        display.setCursor(24,24);
        display.print("Grados");
        display.setCursor(0,32);
        display.print("Pendulo=");
        display.setCursor(0,40);
        display.print("Brazo=");
        display.display();
        pwm_valor=50;
      }
      if(grados_b>270)
      {
        estado=0;
        digitalWrite(motor_b,HIGH);
            //delay(2);
        digitalWrite(motor_a,HIGH);
            
        digitalWrite(ena_motor,HIGH);
        delay(1000);
        detener();
      }else if(grados_b<-270)
      {
        estado=0;
        digitalWrite(motor_b,HIGH);
            //delay(2);
        digitalWrite(motor_a,HIGH);
            
        digitalWrite(ena_motor,HIGH);
        delay(1000);
        detener();
      }
      
      delay(10);
    break;
    
    //Derecha
    case 2:
      if(digitalRead(derecha)==0)
      {
        estado=0;
        digitalWrite(motor_a,HIGH);
        digitalWrite(motor_b,HIGH);
        delay(20);
        digitalWrite(motor_a,LOW);
        digitalWrite(motor_b,LOW);
        digitalWrite(ena_motor,LOW);
        analogWrite(pwm_motor,0);
        digitalWrite(led_derecha,LOW);
      }
      men_grados();
    break;
    
    //Izquierda
    case 3:
      if(digitalRead(izquierda)==0)
      {
        estado=0;
        digitalWrite(motor_a,HIGH);
        digitalWrite(motor_b,HIGH);
        delay(20);
        digitalWrite(motor_a,LOW);
        digitalWrite(motor_b,LOW);
        digitalWrite(ena_motor,LOW);
        analogWrite(pwm_motor,0);
        digitalWrite(led_izquierda,LOW);
      }
      men_grados();
    break;
    
    //Reset
    case 4:
      conteo_b = 32768;
      conteo_p = 32768;
      detener();
      if(digitalRead(reset)==0)
      {
        estado=0;
        delay(50);
      }
    break;    
  }
}

////////////////////////////////
//Interrupciones encoder brazo//
////////////////////////////////

void eb_A() 
{
  if(digitalRead(19)==LOW)
  {
    conteo_b++;
  }else
  {
    conteo_b--;
  }
}

void eb_B() 
{
  if(digitalRead(18)==LOW)
  {
    conteo_b--;
  }else
  {
    conteo_b++;
  }
}

//////////////////////////////////
//Interrupciones encoder pendulo//
//////////////////////////////////

void ep_A() 
{
  if(digitalRead(3)==LOW)
  {
    conteo_p++;
  }else
  {
    conteo_p--;
  }
}

void ep_B() 
{
  if(digitalRead(2)==LOW)
  {
    conteo_p--;
  }else
  {
    conteo_p++;
  }
}

//////////////
//Subrutinas//
//////////////

void men_grados ()
{
  //display.clearDisplay();
  display.fillRect(48, 32, 84, 39, 0);
  //display.display();
  display.setCursor(48,32);
  display.print(grados_p,2);
  //display.display();
  display.fillRect(36, 40, 84, 47, 0);
  //display.display();
  display.setCursor(36,40);
  display.print(grados_b,2);
  display.display();
}

void detener ()
{
  delay(20);
  digitalWrite(motor_a,LOW);
  digitalWrite(motor_b,LOW);
  digitalWrite(ena_motor,LOW);
  pwm_valor=50;
  analogWrite(pwm_motor,0);
  display.fillRect(0, 24, 84, 40, 0);
  display.setCursor(24,24);
  display.print("Grados");
  display.setCursor(0,32);
  display.print("Pendulo=");
  display.setCursor(0,40);
  display.print("Brazo=");
  display.display();
  while (Serial.available())
  {
    String str = Serial.readStringUntil('\n');
    dato = str.toFloat();
  }
}

