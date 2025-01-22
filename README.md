# ACCPPs

Una librería auxiliar de Python para abordar la pregunta **"¿Es posible
predecir péptidos anticancerígenos penetradores de membrana?"** en el **Tanque
de Pensamiento en Ciencia de Datos IIMAS 2025**.

## Instalación

Creamos y activamos un ambiente de conda para instalar ACCPPs y sus dependencias:

```bash
conda create -n accpps python=3.12
conda activate accpps
```

Instalamos las dependencias:

```bash
conda install mamba
mamba install -c conda-forge biopython rdkit
pip install propy3
pip install peptides
```

Clonamos el repositorio de ACCPPs y lo instalamos en el ambiente de conda:
```bash
git clone https://github.com/uibcdf/ACCPPs
cd ACCPPs
pip install --no-deps --editable .
```

## Uso

La librería contiene una función para extraer características físico-químicas
de un péptido dada su secuencia de aminoácidos:

```python
import accpps

sequence = "ACDEFGHIKLMNPQRSTVWY"
features = accpps.get_features(sequence)
```

El objeto de salid `features` es un diccionario con 11097 características físico-químicas.

Además, la librería contiene 4 listas de secuencias de péptidos:

- `accpps.acps`: 1281 secuencias de péptidos anticancerígenos.
- `accpps.non_cps`: 6661 secuencias de péptidos no anticancerígenos.
- `accpps.ccps`: 1582 secuencias de péptidos penetradores de membrana.
- `accpps.non_ccps`: 4369 secuencias de péptidos no penetradores de membrana.

## Advertencia

- Si tu abordaje se basa en el uso de las características físico-químicas, limita
tu base de datos de péptidos a secuencias con más de 4 amino ácidos.

## Consejo

Tener tu propio fork del repositorio e ir implementando tu abordaje como parte
de tu versión de tu librería ACCPPs puede ser muy útil si tu propuesta es exitosa.
