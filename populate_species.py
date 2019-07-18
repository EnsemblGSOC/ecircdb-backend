import requests
from bs4 import BeautifulSoup
from core.models import Species
from django.core.files import File

r = requests.get('http://asia.ensembl.org/info/about/species.html')
soup = BeautifulSoup(r.text, 'html.parser')


def find_species(taxon_id):
    species_list = soup.tbody.find_all('tr')
    for i in range(len(species_list)):
        if species_list[i].find_all('td')[2].text == taxon_id:
            wiki = requests.get(
                'https://en.wikipedia.org/wiki/'+species_list[i].find_all('td')[1].i.text)
            desc_soup = BeautifulSoup(wiki.text, 'html.parser')
            return {
                'thumbnail': 'http://asia.ensembl.org'+species_list[i].td.a.img.get('src'),
                'common_name': species_list[i].td.find_all('a')[1].text,
                'scientific_name': species_list[i].find_all('td')[1].i.text,
                'description': desc_soup.find_all('p')[2].text
            }


for species in Species.objects.all():
    species_details = find_species(str(species.taxon_id))
    species.common_name = species_details['common_name']
    species.scientific_name = species_details['scientific_name']
    species.description = species_details['description']
    thumbnail = requests.get(species_details['thumbnail'])
    if thumbnail.status_code == 200:
        with open("./media/species_thumbnails/"+str(species.taxon_id)+".png", 'wb') as f:
            f.write(thumbnail.content)
        species.thumbnail.save(str(species.taxon_id)+'.png', File(
            open("./media/species_thumbnails/"+str(species.taxon_id)+".png", 'rb')))
    species.save()
