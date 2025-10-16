class House:

	def __init__(self, rprice):
		self._price = rprice

	@property
	def price(self):
	    return self._price
	
	@price.setter
	def price(self, new_price):
            print('set price')
            if new_price > 0 and isinstance(new_price, float):
                self._price = new_price
            else:
                print("Please enter a valid price")

	@price.deleter
	def price(self):
		del self._price

house = House(50000.0)
print(house.price)
