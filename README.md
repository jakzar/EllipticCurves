# Elliptic Curve Cryptography (ECC) in SageMath

Welcome to the Elliptic Curve Cryptography (ECC) repository implemented in SageMath. This repository contains implementations of attacks on DLP on EC and few EC protocols.

## Table of Contents

1. [Introduction](#introduction)
2. [Protocols Included](#protocols-included)
3. [Setup and Requirements](#setup-and-requirements)
4. [How to Run](#how-to-run)
5. [Contributing](#contributing)
6. [License](#license)

## Introduction

//to do :)


## Protocols Included

1. **Custom EC arithmetic test**
    - File: [arithmetic_test.sage](./ec_arithmetic_test.sage)
    - Description: Test of own implementation of arithmetic on EC (mongomery ladder).
  
2. **Rho-method**
    - File: [rho_method.sage](./ec_rho_method.sage)
    - Description: File contains implementation of the rho-method for solving DLP on EC.

3. **Pohlig-Hellman**
    - File: [pohlig_hellman.sage](./ec_pohlig_hellman.sage)
    - Description: Implementation of the Pohlig-Hellman algorithm for solving DLP in subgroups of smaller orders.

4. **ECDSA**
    - File: [ECDSA.sage](./ec_ecdsa.sage)
    - Description: Implementation of the DSA protocol on EC.

5. **ECKCDSA**
    - File: [ECKCDSA.sage](./ec_kcdsa.sage)
    - Description: Implementation of the ECKCDSA protocol.

6. **ECMQV**
    - File: [ECMQV.sage](./ec_ecmqv.sage)
    - Description: Implementation of the ECMQV protocol.


## Setup and Requirements

Before running the protocols, ensure you have SageMath installed on your system. SageMath is available for various platforms, and installation instructions can be found on the official [SageMath website](https://www.sagemath.org/download.html). You can also use [Sage Cell Server](https://sagecell.sagemath.org/)

## How to Run

To run a specific protocol, follow these steps:

1. Open the corresponding SageMath script using your preferred text editor or SageMath environment.
2. Adjust any parameters or configurations if necessary.
3. Execute the script in SageMath.

Feel free to explore and modify the codes to suit your needs...

## Contributing
Contributions to this repository are welcome...

## License
This repository is licensed under the [MIT License](./LICENSE). You are free to use, modify, and distribute the code for both commercial and non-commercial purposes.
