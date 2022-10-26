import hashlib


def hash_32_bit_integer(iterable):
    """Return 32-bit integer hash from an iterable

    :param iterable: iterable containing immutable objects
    :return: hash_32: int, 10-digit integer representing 32-bit hash
    """
    hash = hashlib.sha256()
    for attribute in iterable:
        hash.update(str(attribute).encode())

    hash_32 = int.from_bytes(hash.digest()[:4], byteorder='little')

    return hash_32
