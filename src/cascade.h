class Particle ;
class DatabasePDG2 ;

//class Decayer{
 //DatabasePDG2* database;
 //public:
 //Decayer(DatabasePDG2* _database): database(_database) {};
 //~Decayer() {};
 void decay(Particle *in, int& nprod, Particle** &out);
//};
